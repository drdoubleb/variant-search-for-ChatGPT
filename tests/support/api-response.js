// Minimal adapter for invoking Vercel handlers from tests and Playwright routes.
export function captureResponse() {
    return {
        statusCode: 200,
        headers: {},
        body: null,
        status(code) { this.statusCode = code; return this; },
        setHeader(name, value) { this.headers[name] = value; },
        json(body) { this.body = body; return this; },
        end() { return this; }
    };
}
