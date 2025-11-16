/**
 * BioSynth AI — Vitest/Vite/TS Diagnostic Script
 * Run with:  npx tsx scripts/diagnose.ts  (or npx ts-node scripts/diagnose.ts)
 */

import fs from "fs";
import path from "path";
import { fileURLToPath } from "url";
import stripJsonComments from "strip-json-comments";

const __dirname = path.dirname(fileURLToPath(import.meta.url));
const frontendDir = path.resolve(__dirname, "..");

function exists(p: string) {
  return fs.existsSync(path.resolve(frontendDir, p));
}

function read(p: string) {
  return fs.readFileSync(path.resolve(frontendDir, p), "utf-8");
}

function write(p: string, data: string) {
  fs.writeFileSync(path.resolve(frontendDir, p), data);
}

// ------------------------------
// 1. Verify vitest.config.ts
// ------------------------------
console.log("\n🔍 Checking vitest.config.ts...");

const vitestPath = "vitest.config.ts";
if (!exists(vitestPath)) {
  console.error("❌ vitest.config.ts missing!");
} else {
  let content = read(vitestPath);

  const badImports = [
    "./vite.config",
    "from 'vite'",
    "@vitejs/plugin-react-swc"
  ];

  let fixed = false;
  badImports.forEach(bad => {
    if (content.includes(bad)) {
      console.log(`⚠️  Found problematic import: ${bad}`);
      fixed = true;
    }
  });

  // Ensure plugin-react is used (not SWC)
  if (content.includes("@vitejs/plugin-react-swc")) {
    console.log("⚠️  Replacing SWC plugin with @vitejs/plugin-react...");
    content = content.replace(
      /@vitejs\/plugin-react-swc/g,
      "@vitejs/plugin-react"
    );
    fixed = true;
  }

  // Ensure using vitest/config (not vite)
  if (content.includes("from 'vite'") && content.includes("defineConfig")) {
    console.log("⚠️  Should use 'vitest/config' instead of 'vite'...");
    fixed = true;
  }

  // Ensure plugin typed as any
  if (content.includes("react()") && !content.includes("react() as any")) {
    console.log("⚠️  Plugin should be typed as 'any' to avoid type mismatch...");
    fixed = true;
  }

  if (fixed) {
    console.log("⚠️  vitest.config.ts needs fixes (see warnings above)");
  } else {
    console.log("✅ vitest.config.ts looks correct");
  }
}

// ------------------------------
// 2. Verify test setup file
// ------------------------------
console.log("\n🔍 Checking setup file...");

const setupPath = "src/tests/vitest.setup.ts";
if (!exists(setupPath)) {
  console.log("⚠️  Setup file missing — creating...");
  write(
    setupPath,
    `import "@testing-library/jest-dom";\n`
  );
  console.log("✅ Created setup file");
} else {
  console.log("✅ Setup file OK");
}

// ------------------------------
// 3. Verify tsconfig path aliases
// ------------------------------
console.log("\n🔍 Checking tsconfig.app.json...");

const tsconfigPath = "tsconfig.app.json";
if (exists(tsconfigPath)) {
  // Strip comments from JSON before parsing (TypeScript config allows comments)
  let tsconfigContent = read(tsconfigPath);
  tsconfigContent = stripJsonComments(tsconfigContent);
  const tsconfig = JSON.parse(tsconfigContent);
  let needsUpdate = false;

  tsconfig.compilerOptions = tsconfig.compilerOptions || {};

  // Check types
  const types = tsconfig.compilerOptions.types || [];
  if (!types.includes("vitest/globals") || !types.includes("node")) {
    console.log("⚠️  Adding Vitest types...");
    tsconfig.compilerOptions.types = [
      ...new Set([...types, "vitest/globals", "node"])
    ];
    needsUpdate = true;
  }

  // Check baseUrl
  if (tsconfig.compilerOptions.baseUrl !== ".") {
    console.log("⚠️  Setting baseUrl to '.'...");
    tsconfig.compilerOptions.baseUrl = ".";
    needsUpdate = true;
  }

  // Check paths
  const paths = tsconfig.compilerOptions.paths || {};
  const expectedPaths = {
    "@/*": ["src/*"],
    "@biosynth/engine/*": ["../packages/engine/src/*"]
  };

  let pathsOk = true;
  for (const [key, value] of Object.entries(expectedPaths)) {
    if (JSON.stringify(paths[key]) !== JSON.stringify(value)) {
      pathsOk = false;
      break;
    }
  }

  if (!pathsOk) {
    console.log("⚠️  Updating path aliases...");
    tsconfig.compilerOptions.paths = {
      ...paths,
      ...expectedPaths
    };
    needsUpdate = true;
  }

  if (needsUpdate) {
    write(tsconfigPath, JSON.stringify(tsconfig, null, 2));
    console.log("✅ tsconfig.app.json updated with Vitest globals + aliases");
  } else {
    console.log("✅ tsconfig.app.json looks correct");
  }
} else {
  console.log("❌ tsconfig.app.json missing!");
}

// ------------------------------
// 4. Check required dependencies
// ------------------------------
console.log("\n🔍 Checking dependencies...");

const pkgPath = "package.json";
const pkg = JSON.parse(read(pkgPath));
const deps = { ...pkg.dependencies, ...pkg.devDependencies };

function missing(dep: string) {
  return !deps[dep];
}

const installs: string[] = [];
const uninstalls: string[] = [];

if (missing("@vitejs/plugin-react")) installs.push("@vitejs/plugin-react");
if (missing("vitest")) installs.push("vitest");
if (missing("jsdom")) installs.push("jsdom");
if (deps["@vitejs/plugin-react-swc"]) uninstalls.push("@vitejs/plugin-react-swc");

if (installs.length > 0) {
  console.log("⬆️  Need to install:", installs.join(", "));
}

if (uninstalls.length > 0) {
  console.log("⬇️  Need to remove:", uninstalls.join(", "));
}

if (installs.length === 0 && uninstalls.length === 0) {
  console.log("✅ All dependencies are correct");
}

console.log("\n✅ Diagnostics completed.");

if (installs.length > 0 || uninstalls.length > 0) {
  console.log("👉 Run recommended installs/uninstalls manually:\n");
  if (uninstalls.length > 0) {
    console.log(`npm uninstall ${uninstalls.join(" ")}\n`);
  }
  if (installs.length > 0) {
    console.log(`npm install -D ${installs.join(" ")}\n`);
  }
}

