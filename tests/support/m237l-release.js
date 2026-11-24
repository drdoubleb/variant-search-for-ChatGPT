import assert from 'node:assert/strict';
import { createHash } from 'node:crypto';
import { readFileSync } from 'node:fs';

export const M237L_QUERY = 'NM_000546.6:c.709A>T';
const hash = bytes => createHash('sha256').update(bytes).digest('hex');

export async function checkReleaseAssets(page) {
    for (const [file, selector, attr] of [
        ['script.js', 'script[src*="script.js?v="]', 'src'],
        ['style.css', 'link[href*="style.css?v="]', 'href']
    ]) {
        const expected = hash(readFileSync(new URL(`../../${file}`, import.meta.url)));
        const el = page.locator(selector).first();
        assert.equal(await el.count(), 1, `The page is missing ${file}'s versioned asset reference.`);
        const url = new URL(await el.getAttribute(attr), page.url());
        assert.ok(url.searchParams.get('v')?.endsWith(expected.slice(0, 10)),
            `${file} is referenced at ${url.searchParams.get('v')}. This checkout expects ${expected.slice(0, 10)}.`);
        const res = await page.request.get(url.href, { timeout: 20000 });
        assert.equal(res.status(), 200, `${file} returned HTTP ${res.status()}.`);
        assert.equal(hash(await res.body()), expected, `The bytes served for ${file} differ from this checkout.`);
    }
}

export async function checkM237lRelease(page, timeout = 45000) {
    const tp53 = page.locator('[data-card="tp53-database"]');
    const exact = tp53.getByText('Exact genomic allele matched by assembly, chromosome, position, REF and ALT.', { exact: true });
    await exact.waitFor({ timeout });
    const detail = tp53.getByText('Exact allele details', { exact: true });
    await detail.click();
    const mutationRow = tp53.getByRole('row').filter({ hasText: 'MutationView ID' });
    assert.equal(await mutationRow.getByRole('cell').textContent(), '3082');
    const alleleRow = tp53.getByRole('row').filter({ hasText: 'Genomic allele' });
    assert.match(await alleleRow.getByRole('cell').textContent(), /GRCh37 chr17:7577572 T>A/);
    await detail.click();
    await tp53.getByText('Assay comparison · p.Met237Leu', { exact: true }).waitFor({ timeout });

    const clinvar = page.locator('[data-card="clinvar"]');
    const submissions = clinvar.getByText(/^Submission evidence \(\d+\)$/);
    await submissions.waitFor({ timeout });
    await submissions.click();
    await clinvar.getByRole('heading', { name: /Germline pathogenicity .* Hereditary cancer-predisposing syndrome/ }).waitFor({ timeout });
    await clinvar.getByText(/Color Health.*Counts toward ClinVar aggregate/).waitFor({ timeout });
    await submissions.click();
}
