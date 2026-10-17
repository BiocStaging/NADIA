// Runs the Excel export of inst/js/nadia-tables/nadia-tables.js outside a
// browser, with ExcelJS, PapaParse and reactable replaced by small stand-ins,
// and prints the values written to the data rows as JSON.
//
// node excel-harness.js <nadia-tables.js> <settings.json> <table.tsv>
const fs = require("fs");
const vm = require("vm");
const [script, settingsFile, tsvFile] = process.argv.slice(2);

const cells = {};
const cell = (r, c) => (cells[r + "," + c] = cells[r + "," + c] || {});
const sheet = {
  getColumn: () => ({}),
  getCell: cell,
  mergeCells: () => {},
  getRow: (r) => ({ getCell: (c) => cell(r, c), commit: () => {} }),
  views: [],
};
const tsv = fs.readFileSync(tsvFile, "utf8");
const ctx = {
  console,
  Date,
  window: { URL: { createObjectURL: () => "blob:", revokeObjectURL: () => {} } },
  document: { createElement: () => ({ click: () => {} }) },
  Blob: function () {},
  Reactable: { getDataCSV: () => tsv },
  // The test table has no quoted fields, so a plain split is enough
  Papa: {
    parse: (text) => {
      const lines = text.trim().split("\n").map((l) => l.split("\t"));
      const fields = lines[0];
      const data = lines.slice(1).map((v) => Object.fromEntries(fields.map((f, i) => [f, v[i]])));
      return { data, meta: { fields } };
    },
  },
  ExcelJS: {
    Workbook: function () {
      this.addWorksheet = () => sheet;
      this.xlsx = { writeBuffer: async () => new Uint8Array(0) };
    },
  },
};
vm.createContext(ctx);
vm.runInContext(fs.readFileSync(script, "utf8"), ctx);

ctx.window.NadiaTables.register("t", JSON.parse(fs.readFileSync(settingsFile, "utf8")));
ctx.window.NadiaTables.exportExcel("t").then(() => {
  const rows = {};
  for (const [key, value] of Object.entries(cells)) {
    const [r, c] = key.split(",").map(Number);
    if (r >= 3 && "value" in value) (rows[r] = rows[r] || {})[c] = value.value;
  }
  console.log(JSON.stringify(rows));
});
