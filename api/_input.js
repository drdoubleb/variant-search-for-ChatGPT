// Strict, low-cost validation for public serverless query parameters.
// Vercel req.query values are strings, but repeated keys can produce arrays.
// Reject those instead of silently coercing attacker-controlled objects.
export function boundedQueryText(value, maxLength = 256) {
    if (typeof value !== 'string') return null;
    const text = value.trim();
    return text.length <= maxLength ? text : null;
}

export function positiveQueryInteger(value, max = 1_000_000_000) {
    const text = boundedQueryText(value, 12);
    if (text === null || !/^[1-9][0-9]*$/.test(text)) return null;
    const number = Number(text);
    return Number.isSafeInteger(number) && number <= max ? number : null;
}

export function validatedChromosome(value) {
    const text = boundedQueryText(value, 6);
    if (!text) return null;
    const chrom = text.replace(/^chr/i, '').toUpperCase();
    return /^(?:[1-9]|1[0-9]|2[0-2]|X|Y|M|MT)$/.test(chrom) ? chrom : null;
}

export function validatedGene(value) {
    const text = boundedQueryText(value, 64);
    return text && /^[A-Za-z0-9_.-]+$/.test(text) ? text : null;
}

export function validatedAllele(value, maxLength = 1000) {
    const text = boundedQueryText(value, maxLength);
    return text && /^[ACGTN]+$/i.test(text) ? text.toUpperCase() : null;
}
