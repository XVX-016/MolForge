export const BASE_URL = import.meta.env.VITE_BACKEND_URL || "http://127.0.0.1:8000";

export const api = {
    predict: (payload: any) =>
        fetch(`${BASE_URL}/api/predict/property`, {
            method: "POST",
            headers: { "Content-Type": "application/json" },
            body: JSON.stringify(payload),
        }).then(r => r.json()),

    predictBatch: (payload: any) =>
        fetch(`${BASE_URL}/api/predict/batch`, {
            method: "POST",
            headers: { "Content-Type": "application/json" },
            body: JSON.stringify(payload),
        }).then(r => r.json()),

    relax: (molfile: string) =>
        fetch(`${BASE_URL}/api/molecule/relax`, {
            method: "POST",
            headers: { "Content-Type": "application/json" },
            body: JSON.stringify({ molfile }),
        }).then(r => r.json()),

    generate3D: (smiles: string) =>
        fetch(`${BASE_URL}/api/molecule/generate-3d`, {
            method: "POST",
            headers: { "Content-Type": "application/json" },
            body: JSON.stringify({ smiles }),
        }).then(r => r.json()),

    train: () =>
        fetch(`${BASE_URL}/api/models/train`, { method: "POST" }),

    listModels: () =>
        fetch(`${BASE_URL}/api/models/list`).then(r => r.json()),
};
