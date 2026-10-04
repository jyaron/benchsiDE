import os, re, json, subprocess, time
JSC = "/System/Library/Frameworks/JavaScriptCore.framework/Versions/Current/Helpers/jsc"
STUB = r"""
var plotlyCalls = [];
function makeEl(id){
  const t = {id:id, style:{}, dataset:{}, value:"", checked:false, _c:"", open:false, disabled:false,
    get innerHTML(){ return this._c; }, set innerHTML(v){ this._c=String(v); },
    get textContent(){ return this._c.replace(/<[^>]+>/g,""); }, set textContent(v){ this._c=String(v).replace(/&/g,"&amp;").replace(/</g,"&lt;"); },
    clientWidth:620, clientHeight:340, children:[], classList:{add(){},remove(){},toggle(){},contains(){return false;}},
    _h:{}, options:[], add(o){ this.options.push(o); }, addEventListener(ev,fn){ this._h[ev]=fn; }, removeEventListener(){}, appendChild(c){ this.children.push(c); return c; }, remove(){},
    insertAdjacentHTML(){}, setAttribute(){}, getAttribute(){return null;}, querySelector(){return makeEl("_q");}, querySelectorAll(){return [];},
    click(){}, focus(){}, scrollIntoView(){}, getContext(){ return {measureText:s=>({width:String(s).length*7}), font:""}; },
    getBoundingClientRect(){ return {width:620,height:340,top:0,left:0}; }, on(ev,f){ (this._on=this._on||{})[ev]=f; }, removeAllListeners(){} };
  return t;
}
var elCache = {};
var document = { getElementById(id){ if(!elCache[id]) elCache[id]=makeEl(id); return elCache[id]; },
  querySelector(){ return makeEl("_q"); }, querySelectorAll(){ return []; }, createElement(t){ return makeEl("_c_"+t); },
  addEventListener(){}, body:makeEl("body"), documentElement:makeEl("html") };
var window = { addEventListener(){}, innerWidth:1400, matchMedia(){ return {matches:false, addEventListener(){}}; }, location:{search:"",hash:""} };
var navigator = { userAgent:"jsc", clipboard:{writeText(){ return Promise.resolve(); }} };
var localStorage = { getItem(){return null;}, setItem(){}, removeItem(){} };
var URL = { createObjectURL(){ return "blob:"; }, revokeObjectURL(){} };
var LASTBLOB=null; function Blob(p,o){ this.p=p; this.o=o; LASTBLOB=this; }
var alert=function(m){ print("ALERT: "+m); }, confirm=function(){return true;};
var __timers=[]; var setTimeout = function(fn){ __timers.push(fn); return __timers.length; }, clearTimeout = function(){};
function __drain(max){ let n=0; for(;;){ if(typeof drainMicrotasks==="function") drainMicrotasks(); if(!__timers.length||n>=(max||1e7)) break; const f=__timers.shift(); f(); n++; } if(typeof drainMicrotasks==="function") drainMicrotasks(); return n; }
var requestAnimationFrame = function(fn){ __timers.push(fn); };
var console = { log(){}, warn(){}, error(){} };
var performance = { now(){ return Date.now(); } };
function Option(t,v){ this.text=t; this.value=v; this.textContent=t; }
var HTMLElement=function(){};
var Plotly = { react(div,tr,lay,cfg){ plotlyCalls.push({id:(div&&div.id)||div, traces:tr, layout:lay, cfg:cfg}); return Promise.resolve(); },
  newPlot(div,tr,lay,cfg){ return Plotly.react(div,tr,lay,cfg); }, purge(){}, relayout(){ return Promise.resolve(); },
  toImage(){ return Promise.resolve("data:"); }, downloadImage(){ return Promise.resolve(); }, Plots:{resize(){}} };
"""
def app_of(src):
    return re.findall(r"<script>\n(.*?)\n</script>", src, re.S)[-1]
def run(app_src, body, extra="", timeout=900):
    """Run app_src + extra + body in JavaScriptCore with the DOM/Plotly stub. Uses a unique temp file (parallel-safe)."""
    import tempfile
    fd, path = tempfile.mkstemp(suffix=".js", dir=os.getcwd()); os.close(fd)
    open(path,"w").write(STUB + app_src + "\n" + extra + "\n__drain();\ntry {\n" + body + "\n__drain();\n} catch(e) { print('RUNTIME FAIL: ' + e + '\\n' + (e.stack||'')); }\n")
    t0 = time.time()
    try:
        o = subprocess.run([JSC, path], capture_output=True, text=True, timeout=timeout)
    finally:
        os.remove(path)
    return o.stdout + o.stderr + f"\n[wall {time.time()-t0:.1f}s]"

def lines_with(out, prefix):
    """Return JSON payloads of stdout lines starting with prefix (e.g. 'JSON ')."""
    return [json.loads(l[len(prefix):]) for l in out.split("\n") if l.startswith(prefix)]

# Typical session setup (count data, human): 
# el("species").value="human"; useBuiltinAnno(); el("normsel").value="tmm"; el("fmode").value="fbe";
# loadFromText(MATRIX_TEXT, DESIGN_TEXT); applyFactor("<factor>"); GROUPTYPE="auto"; analyzeNow(); __drain();
# el("selA").value="<ref>"; el("selB").value="<test>"; el("selFDR").value="0.05"; el("fcthr").value="1";
# deMethod="mod"|"voom"|"welch"; DECOVSEL=[...covariate names]; deKey=""; const r=computeDE();  // r.lfc, r.tv, r.p, r.q, r.se, r.aM, r.bM, r.prior
# Globals after analyzeNow: nG, nS, GENES, IDS, L (log2 matrix, gene-major, L[g*nS+q]), S (samples), included, CNT, TMMF, EFFLIB, GROUPS, groupsIdx()
