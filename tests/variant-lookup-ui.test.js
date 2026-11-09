// Exercise the real form handler and rendered cards with live-response fixtures.
// LinkeDOM supplies the DOM; every external request is mocked (no live services).
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import vm from 'node:vm';
import { parseHTML } from 'linkedom';
const script = readFileSync(new URL('../script.js', import.meta.url), 'utf8');
const html = readFileSync(new URL('../index.html', import.meta.url), 'utf8');
const fixture = JSON.parse(readFileSync(new URL('./fixtures/pten-utr.json', import.meta.url)));

const brafAnnotation = {
    _id:'chr7:g.140453136A>T', dbnsfp:{genename:'BRAF'},
    snpeff:{ann:{feature_id:'NM_004333.4',hgvs_c:'c.1799T>A',hgvs_p:'p.Val600Glu',genename:'BRAF',effect:'missense_variant'}},
    clinvar:{variant_id:13961,gene:{symbol:'BRAF'},rcv:{clinical_significance:'Pathogenic'}}
};
async function lookup(input, { vepDown = false, maneDown = false, extended = false, braf = false } = {}) {
    const { document } = parseHTML(html);
    const form = document.getElementById('variantForm');
    const listen = form.addEventListener.bind(form);
    let submit;
    form.addEventListener = (event, fn, ...args) => {
        if (event === 'submit') submit = fn;
        return listen(event, fn, ...args);
    };
    const requests = [], errors = [], timers = [];
    const response = (data, status = 200) => ({ ok: status === 200, status,
        json: async () => structuredClone(data), text: async () => JSON.stringify(data) });
    const fetch = async address => {
        const url = new URL(address);
        requests.push(url);
        if (braf && url.hostname === 'myvariant.info') {
            if (url.pathname.startsWith('/v1/variant/')) {
                return decodeURIComponent(url.pathname).endsWith(brafAnnotation._id)
                    ? response(brafAnnotation) : response({},404);
            }
            return response({hits:[brafAnnotation]});
        }
        if (braf && url.pathname.includes('/variant_recoder/')) return response({},400);
        if (braf && url.pathname === '/api/clinvar-variant') return response({germline:{description:'Pathogenic'}});
        if (braf && url.pathname === '/api/clinvar-region') return response({variants:[],total:0});
        if (url.hostname === 'myvariant.info') {
            if (url.pathname.startsWith('/v1/variant/')) return response({ error:'not found' },404);
            const q = url.searchParams.get('q');
            const key = q === 'PTEN c.-8T>G' ? 'mv-gene-space'
                : q === 'NM_000314.8 c.-8T>G' ? 'mv-nm-short-space'
                : q === 'NM_001304717.5 c.512T>G' ? 'mv-nm-long-space' : null;
            return response({ hits: key ? fixture.searchHits[key] : [] });
        }
        if (url.pathname.includes('/variant_recoder/')) {
            const result = structuredClone(fixture.recoder);
            result[0].G.input = decodeURIComponent(url.pathname.split('/').at(-1));
            return response(result);
        }
        if (url.pathname.includes('/lookup/symbol/')) return maneDown || braf ? response({},404) : response(fixture.mane);
        if (url.pathname.includes('/vep/')) return vepDown ? response({},400) : response(fixture.vep);
        if (url.pathname === '/api/clinvar-region') {
            assert.equal(url.searchParams.get('pos'), '89624219', 'region must retain original locus');
            return response({ variants: [fixture.clinvar], total:1 });
        }
        if (url.pathname === '/api/clinvar-variant') {
            assert.equal(url.searchParams.get('id'), '4757098', 'never request unrelated ClinVar IDs');
            return response({ variationId:'4757098', germline:{description:'Uncertain significance'} });
        }
        // Unrelated cards fail gracefully; they must not affect variant identity.
        return response({},404);
    };
    const window = { location: new URL('https://drdoubleb.com/variant/'),
        history:{ replaceState() {} }, matchMedia:()=>({matches:false}) };
    const context = vm.createContext({ document, window, fetch,
        console: { log(){}, warn(){}, error(e){errors.push(String(e));} },
        URL, URLSearchParams, AbortController, performance,
        setTimeout(fn, ms) { const timer = setTimeout(fn, Math.min(ms, 5)); timers.push(timer); return timer; },
        clearTimeout
    });
    try {
        vm.runInContext(script, context);
        vm.runInContext('initVariantSearchUi()', context);
        document.getElementById('variantInput').value = input;
        await submit({preventDefault(){}});
        await new Promise(resolve => setTimeout(resolve, 15));
        assert.equal(document.getElementById('status').textContent, '', errors.join('\n'));
        assert.equal(document.getElementById('resultSection').classList.contains('hidden'),false);
        const card = name => [...document.querySelectorAll('.card')].find(el=>el.querySelector('h3')?.textContent===name);
        const variant = card('Variant');
        assert(variant, 'variant card rendered');
        const expanded = variant.querySelector('details');
        const transcriptDetails = expanded?.textContent || '';
        expanded?.remove();
        const main = variant.textContent;
        if (braf) {
            assert.match(main,/chr7:g\.140453136A>T/);
            assert.match(main,/7-140453136-A-T/);
            assert.match(main,/Transcript: NM_004333\.4/);
            assert.match(main,/c\.: c\.1799T>A/);
            assert.match(main,/p\.: p\.Val600Glu/);
            assert.match(card('ClinVar').textContent,/Variation ID: 13961/);
            console.log(`ok UI lookup: ${input} (recoder unavailable)`);
            return;
        }
        assert.match(main, /chr10:g\.89624219T>G/);
        assert.match(main, /10-89624219-T-G/);
        if (extended) {
            assert.match(main,/Transcript: NM_001304717\.5/);
            assert.match(main,/c\.: c\.512T>G/);
            assert.match(main,/p\.: p\.Leu171Arg/);
            assert.match(main,/missense_variant/);
        } else {
            assert.match(main, /Transcript: NM_000314\.8/);
            if (!maneDown) assert.match(main,/\(MANE Select\)/);
            assert.match(main, /c\.: c\.-8T>G/);
            assert.match(main, /p\.: N\/A/);
            assert.match(main, /5_prime_UTR_variant/);
            assert.doesNotMatch(main, /Leu171Arg|NON_SYNONYMOUS|897119/);
            if (!vepDown) assert.match(transcriptDetails, /NM_001304717\.5: c\.512T>G, NP_001291646\.4:p\.Leu171Arg/);
            else assert.doesNotMatch(transcriptDetails, /Leu171Arg/);
        }
        const clinvar = card('ClinVar');
        assert.match(clinvar.textContent,/Variation ID: 4757098/);
        assert.match(clinvar.textContent,/Uncertain significance/);
        assert(clinvar.querySelector('a[href="https://www.ncbi.nlm.nih.gov/clinvar/variation/4757098/"]'));
        assert(!requests.some(url=>url.pathname.includes('/map/')), 'hg19 misses must not trigger speculative liftover');
        if (!extended) assert(!requests.some(url=>url.pathname==='/api/clinvar-protein'), 'UTR view must not search another isoform\'s protein');
        console.log(`ok UI lookup: ${input}${vepDown ? ' (VEP unavailable)' : ''}`);
    } finally {
        timers.forEach(clearTimeout);
    }
}
await lookup('PTEN c.-8T>G');
await lookup('10 89624219 PTEN T G');
await lookup('NM_000314.8:c.-8T>G');
await lookup('10 89624219 PTEN T G', {vepDown:true});

await lookup('NM_001304717.5:c.512T>G', {extended:true});
await lookup('10 89624219 PTEN T G', {maneDown:true});
await lookup('BRAF V600E', {braf:true});
await lookup('7 140453136 BRAF A T', {braf:true});
await lookup('NM_004333.4:c.1799T>A', {braf:true});
