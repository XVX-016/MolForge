import axios from "axios";

const BASE_URL = "http://localhost:8000"; // Force local for debugging

export const apiClient = axios.create({
    baseURL: BASE_URL,
    headers: {
        "Content-Type": "application/json",
    },
});
