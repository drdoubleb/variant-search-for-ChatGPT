import { SaxesParser } from 'saxes';

// EFetch XML has mixed text/markup. Keep document order so decoding inline
// elements in an abstract or submission comment does not change its wording.
export function parseXml(xml, { allowPubMedDoctype = false } = {}) {
    const document = { children: [], content: [] };
    const stack = [document];
    const parser = new SaxesParser();
    parser.on('doctype', declaration => {
        // PubMed emits a reference to its public DTD. Saxes never resolves it.
        // Deny internal subsets and all other declarations.
        if (!allowPubMedDoctype || !/^\s*PubmedArticleSet\s+(?:PUBLIC|SYSTEM)\s+[^\[\]]+$/.test(declaration)) {
            throw new Error('Unexpected XML document type');
        }
    });
    parser.on('error', () => { throw new Error('Malformed upstream XML'); });
    parser.on('opentag', ({ name, attributes }) => {
        const parent = stack.at(-1);
        const element = { name, attributes, children: [], content: [] };
        parent.children.push(element);
        parent.content.push(element);
        stack.push(element);
    });
    parser.on('closetag', () => stack.pop());
    const appendText = value => stack.at(-1).content.push(value);
    parser.on('text', appendText);
    parser.on('cdata', appendText);
    parser.write(xml).close();
    return document.children[0];
}

export const children = (node, name) => node?.children.filter(element => element.name === name) ?? [];
export const child = (node, name) => node?.children.find(element => element.name === name);
export const attribute = (node, name) => node?.attributes[name] ?? '';

function innerText(node) {
    return node.content.map(part => typeof part === 'string' ? part : innerText(part)).join('');
}

export const text = node => node ? innerText(node).trim() : '';
export const childText = (node, name) => text(child(node, name));
