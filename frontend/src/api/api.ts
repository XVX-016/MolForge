import axios from "axios";

const BASE_URL = "http://127.0.0.1:8000"; // Use 127.0.0.1 instead of localhost for Windows reliability

export const apiClient = axios.create({
    baseURL: BASE_URL,
    headers: {
        "Content-Type": "application/json",
    },
});
