"""End-to-end checks for Cactus Alignment Viewer.

1. python tests/make_test_data.py     (writes tests/data/: bigMaf, MAF, bigWig, bigBed, FASTA+fai, GFF3, tree)
2. python tests/run_tests.py          (opens index.html in headless Chromium and checks the readers)

Requires: pip install playwright pybigtools pyBigWig && python -m playwright install chromium
"""
import json, pathlib, sys, pyBigWig
from playwright.sync_api import sync_playwright

ROOT = pathlib.Path(__file__).resolve().parent.parent
DATA = ROOT / "tests" / "data"
APP = (ROOT / "index.html").as_uri()
fails = []

def check(name, ok, detail=""):
    print(("PASS " if ok else "FAIL ") + name + (f"  {detail}" if detail else ""))
    if not ok: fails.append(name)

truth = json.load(open(DATA / "truth.json"))
ref = "".join(open(DATA / "ref.fa").read().split(">")[1].split("\n")[1:])
with sync_playwright() as p:
    b = p.chromium.launch(); pg = b.new_page(viewport={"width": 1500, "height": 950}, accept_downloads=True)
    errors = []; pg.on("pageerror", lambda e: errors.append(str(e)))
    pg.goto(APP); pg.wait_for_timeout(300)
    files = ["test.bigMaf.bb", "tree.nwk", "atac.bw", "peaks.bb", "plain.bed", "genes.gff3", "ref.fa", "ref.fa.fai"]
    pg.set_input_files("#fileInput", [str(DATA / f) for f in files])
    pg.wait_for_function("S.tracks.length>=6", timeout=60000); pg.wait_for_timeout(1500)
    res = pg.evaluate("""async (truth)=>{const t=mafTrack();let ok=0,bad=0;const by={};
      for(const k of Object.keys(truth)){const [sp,pos]=k.split('|');(by[Math.floor(+pos/2000)]=by[Math.floor(+pos/2000)]||[]).push([sp,+pos,truth[k]]);}
      for(const w of Object.keys(by)){const s=+w*2000;const d=await computeMafBins(t,'NC_069100.1',s,s+2000,2000,{noCache:true});
        for(const [sp,pos,ch] of by[w]){const a=d.sp.get(sp);const c=a&&a.chars[pos-s];const got=c?String.fromCharCode(c):(a&&a.gap[pos-s]?'-':'?');if(got===ch)ok++;else bad++;}}
      return{ok,bad};}""", truth)
    check("bigMaf bases match ground truth", res["bad"] == 0 and res["ok"] > 1000, str(res))
    bw = pyBigWig.open(str(DATA / "atac.bw"))
    got = pg.evaluate("(async()=>{const t=S.tracks.find(x=>x.kind==='signal');return [...(await t.bbi.wigBins('NC_069100.1',100000,100100,4)).max]})()")
    exp = [bw.stats("NC_069100.1", 100000 + i * 25, 100025 + i * 25, type="max")[0] for i in range(4)]
    check("bigWig values match pyBigWig", all(abs(a - b) < 1e-4 for a, b in zip(got, exp)))
    fa = pg.evaluate("(async()=>{const A=await getAlignmentColumns({chrom:'NC_069100.1',start:150000,end:150300},{insertions:false});return A.seqs[A.ids.indexOf(S.refId)]})()")
    check("exported reference row equals genome FASTA", fa.upper() == ref[150000:150300].upper())
    pg.set_input_files("#fileInput", [str(DATA / "test.maf.gz")]); pg.wait_for_function("mafTrack().aln.constructor.name==='TextMafSource'", timeout=60000)
    res2 = pg.evaluate("""async (truth)=>{const t=mafTrack();let ok=0,bad=0;for(const k of Object.keys(truth).slice(0,400)){const [sp,pos]=k.split('|');const d=await computeMafBins(t,'NC_069100.1',+pos,+pos+1,1);
      const a=d.sp.get(sp);const got=a?(a.chars[0]?String.fromCharCode(a.chars[0]):a.gap[0]?'-':'?'):'?';if(got===truth[k])ok++;else bad++;}return{ok,bad};}""", truth)
    check("text MAF (.maf.gz) bases match ground truth", res2["bad"] == 0, str(res2))
    check("no JavaScript errors", not errors, "; ".join(errors[:3]))
    b.close()
sys.exit(1 if fails else 0)
