/*
 * Tests for three-letter -> single-letter protein HGVS conversion (helpers copied
 * verbatim from script.js — KEEP IN SYNC). Beyond simple substitutions, the summary
 * card must also render single-letter codes for in-frame deletions, duplications,
 * insertions and delins so that e.g. PARP3 p.Lys18del is shown as p.K18del. This
 * matters for the Ensembl VEP fallback, which returns three-letter RefSeq protein
 * notations (NP_005476.4:p.Lys18del) when MyVariant.info has no annotation.
 *
 * Run with: node tests/protein-notation.test.js
 */

// --- helpers copied verbatim from script.js -------------------------------

function tripleToSingle(prot) {
    if (!prot) return null;
    const m = prot.match(/^([A-Z]{3})(\d+)([A-Z]{3})$/);
    if (!m) return null;
    const aaSingle = {
        ALA: 'A', ARG: 'R', ASN: 'N', ASP: 'D', CYS: 'C', GLN: 'Q', GLU: 'E', GLY: 'G',
        HIS: 'H', ILE: 'I', LEU: 'L', LYS: 'K', MET: 'M', PHE: 'F', PRO: 'P', SER: 'S',
        THR: 'T', TRP: 'W', TYR: 'Y', VAL: 'V'
    };
    const ref = aaSingle[m[1]];
    const pos = m[2];
    const alt = aaSingle[m[3]];
    if (ref && alt) return `${ref}${pos}${alt}`;
    return null;
}

function aaThreeToSingle(aa) {
    if (!aa) return null;
    const aaSingle = {
        ALA: 'A', ARG: 'R', ASN: 'N', ASP: 'D', CYS: 'C', GLN: 'Q', GLU: 'E', GLY: 'G',
        HIS: 'H', ILE: 'I', LEU: 'L', LYS: 'K', MET: 'M', PHE: 'F', PRO: 'P', SER: 'S',
        THR: 'T', TRP: 'W', TYR: 'Y', VAL: 'V',
        TER: '*', STOP: '*'
    };
    return aaSingle[String(aa).toUpperCase()] || null;
}

function convertProteinBodyToSingle(proteinBody) {
    if (!proteinBody) return null;
    const body = String(proteinBody).trim();

    // Frameshift with explicit downstream stop: Val133GlyfsTer47 -> V133Gfs*47
    let m = body.match(/^([A-Za-z]{3})(\d+)([A-Za-z]{3})fs(?:Ter|Stop|\*)(\d+)$/i);
    if (m) {
        const ref = aaThreeToSingle(m[1]);
        const alt = aaThreeToSingle(m[3]);
        if (ref && alt) return `${ref}${m[2]}${alt}fs*${m[4]}`;
    }

    // Frameshift without explicit alt amino acid: Arg97fsTer12 -> R97fs*12
    m = body.match(/^([A-Za-z]{3})(\d+)fs(?:Ter|Stop|\*)(\d+)$/i);
    if (m) {
        const ref = aaThreeToSingle(m[1]);
        if (ref) return `${ref}${m[2]}fs*${m[3]}`;
    }

    // Simple substitution / nonsense: Arg714Ter -> R714*
    m = body.match(/^([A-Za-z]{3})(\d+)([A-Za-z]{3}|Ter|Stop|\*)$/i);
    if (m) {
        const ref = aaThreeToSingle(m[1]);
        const alt = aaThreeToSingle(m[3]) || (m[3] === '*' ? '*' : null);
        if (ref && alt) return `${ref}${m[2]}${alt}`;
    }

    // In-frame deletions, duplications, insertions and delins (single residue or range):
    //   Lys18del -> K18del, Lys18_Arg20del -> K18_R20del,
    //   Lys18_Arg20delinsGlyPro -> K18_R20delinsGP, Met1? -> unchanged.
    if (/(?:del|dup|ins)/i.test(body) && /[A-Za-z]{3}\d/.test(body)) {
        const converted = body.replace(
            /Ala|Arg|Asn|Asp|Cys|Gln|Glu|Gly|His|Ile|Leu|Lys|Met|Phe|Pro|Ser|Thr|Trp|Tyr|Val|Ter|Stop/gi,
            (token) => aaThreeToSingle(token) || token
        );
        if (converted !== body) return converted;
    }

    // Fallback to strict triple substitution converter.
    const compact = body.toUpperCase();
    return tripleToSingle(compact);
}

function formatProteinDisplayWithSingleLetter(proteinHgvs) {
    if (!proteinHgvs) return '';
    const proteinText = String(proteinHgvs).trim();
    const proteinIndex = proteinText.toLowerCase().lastIndexOf('p.');
    if (proteinIndex === -1) return proteinText;
    const proteinBody = proteinText.slice(proteinIndex + 2).trim();
    const single = convertProteinBodyToSingle(proteinBody);
    if (!single) return proteinText;
    return `${proteinText} (p.${single})`;
}

// --- assertions -----------------------------------------------------------

let passed = 0;
let failed = 0;
function check(name, cond) {
    if (cond) {
        passed++;
    } else {
        failed++;
        console.error(`FAIL: ${name}`);
    }
}

// Substitutions and nonsense (pre-existing behaviour must be preserved).
check('substitution Val600Glu -> V600E', convertProteinBodyToSingle('Val600Glu') === 'V600E');
check('substitution Arg273His -> R273H', convertProteinBodyToSingle('Arg273His') === 'R273H');
check('nonsense Arg714Ter -> R714*', convertProteinBodyToSingle('Arg714Ter') === 'R714*');
check('frameshift Val133GlyfsTer47 -> V133Gfs*47', convertProteinBodyToSingle('Val133GlyfsTer47') === 'V133Gfs*47');

// In-frame deletions / duplications / insertions / delins (new behaviour).
check('single-residue del Lys18del -> K18del', convertProteinBodyToSingle('Lys18del') === 'K18del');
check('range del Lys18_Arg20del -> K18_R20del', convertProteinBodyToSingle('Lys18_Arg20del') === 'K18_R20del');
check('single-residue dup Lys18dup -> K18dup', convertProteinBodyToSingle('Lys18dup') === 'K18dup');
check('range dup Ala767_Val769dup -> A767_V769dup', convertProteinBodyToSingle('Ala767_Val769dup') === 'A767_V769dup');
check('delins Gly12delinsAlaSer -> G12delinsAS', convertProteinBodyToSingle('Gly12delinsAlaSer') === 'G12delinsAS');
check('range delins Lys18_Arg20delinsGlyPro -> K18_R20delinsGP',
    convertProteinBodyToSingle('Lys18_Arg20delinsGlyPro') === 'K18_R20delinsGP');

// End-to-end display: the summary card appends the single-letter form in parentheses.
check('display p.Lys18del -> p.Lys18del (p.K18del)',
    formatProteinDisplayWithSingleLetter('p.Lys18del') === 'p.Lys18del (p.K18del)');
check('display p.Val600Glu -> p.Val600Glu (p.V600E)',
    formatProteinDisplayWithSingleLetter('p.Val600Glu') === 'p.Val600Glu (p.V600E)');

console.log(`\nProtein-notation tests: ${passed} passed, ${failed} failed`);
if (failed > 0) process.exit(1);
