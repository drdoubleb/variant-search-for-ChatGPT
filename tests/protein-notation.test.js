/*
 * Tests for three-letter -> single-letter protein HGVS conversion (helpers imported
 * from script.js). Beyond simple substitutions, the summary
 * card must also render single-letter codes for in-frame deletions, duplications,
 * insertions and delins so that e.g. PARP3 p.Lys18del is shown as p.K18del. This
 * matters for the Ensembl VEP fallback, which returns three-letter RefSeq protein
 * notations (NP_005476.4:p.Lys18del) when MyVariant.info has no annotation.
 *
 * Run with: node tests/protein-notation.test.js
 */

// --- real helpers imported from script.js ---------------------------------

await import('../script.js');
const { convertProteinBodyToSingle, formatProteinDisplayWithSingleLetter } = globalThis.__variantSearchHelpers;

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
