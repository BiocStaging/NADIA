/*
 * NADIA results tables: filter panel, condition chips and Excel export.
 *
 * Shared by results_list_widget(), protein_list_widget(), quant_list_widget()
 * and summary_list_widget(). Each widget wraps its controls and its table in an
 * element carrying data-nadia-table="<elementId>" and registers its settings
 * with NadiaTables.register(); every function then takes the elementId, so two
 * tables on one page never act on each other.
 *
 * Needs reactable (Reactable.*), and for the export ExcelJS and PapaParse,
 * which the package ships alongside this file.
 */
(function () {
  'use strict';

  var tables = {};

  function root(id) {
    var nodes = document.querySelectorAll('[data-nadia-table]');
    for (var i = 0; i < nodes.length; i++) {
      if (nodes[i].getAttribute('data-nadia-table') === id) return nodes[i];
    }
    return document;
  }

  function chips(id) {
    return root(id).querySelectorAll('.pl-cond-chip');
  }

  function hex2argb(hex) {
    var h = (hex || '').replace('#', '');
    if (h.length !== 6) return 'FF808080';
    return 'FF' + h.toUpperCase();
  }

  // A number only when the whole value is one: Number("60S ribosomal protein")
  // is NaN, where parseFloat() would keep the 60
  function toNumber(v) {
    var s = String(v).trim();
    if (s === '') return null;
    var n = Number(s);
    return isFinite(n) ? n : null;
  }

  function solid(argb) {
    return { type: 'pattern', pattern: 'solid', fgColor: { argb: argb } };
  }

  function thinBorder(argb) {
    var side = argb ? { style: 'thin', color: { argb: argb } } : { style: 'thin' };
    return { top: side, left: side, bottom: side, right: side };
  }

  // Shows or hides the conditions, as columns or as rows
  function applyConditions(id) {
    var t = tables[id];
    var hidden = Object.keys(t.hidden).filter(function (c) { return t.hidden[c]; });
    if (t.conditionMode === 'rows') {
      Reactable.setFilter(id, 'Condition', hidden.length > 0 ? hidden.join(',') : undefined);
    } else {
      var cols = [];
      hidden.forEach(function (c) {
        if (t.conditions[c]) cols = cols.concat(t.conditions[c]);
      });
      Reactable.setHiddenColumns(id, cols);
    }
  }

  // --- Excel layouts --------------------------------------------------------

  // One header row, renamed, over the columns the table exports (DE results)
  function sheetFlat(ws, rows, fields, x) {
    ws.columns = fields.map(function (h) {
      return { header: x.headers[h] || h, key: h, width: x.widths[h] || x.defaultWidth };
    });
    var headerRow = ws.getRow(1);
    headerRow.font = { bold: true, color: { argb: 'FFFFFFFF' }, size: 11 };
    headerRow.fill = solid(x.headerFill);
    headerRow.alignment = { vertical: 'middle', horizontal: 'center' };
    headerRow.height = 30;

    rows.forEach(function (row) {
      var rowData = {};
      fields.forEach(function (h) { rowData[h] = row[h]; });
      var addedRow = ws.addRow(rowData);
      x.numeric.forEach(function (col) {
        if (fields.indexOf(col) === -1) return;
        var cell = addedRow.getCell(col);
        if (cell && cell.value) {
          var num = parseFloat(cell.value);
          if (!isNaN(num)) cell.value = num;
        }
      });
    });

    ws.eachRow(function (row) {
      row.eachCell(function (cell) { cell.border = thinBorder(); });
    });
  }

  // Merged group headers over per-sample headers tinted by condition
  // (Protein_ID and Protein_QUANT). The annotation columns listed in x.text
  // (identifiers, descriptions, genes) are written as text, whatever they
  // start with.
  function sheetGrouped(ws, rows, x) {
    var text = x.text || [];
    var flatCols = [];
    var flatHeaders = [];
    var colCondMap = {};
    x.groups.forEach(function (g) {
      g.columns.forEach(function (c, i) {
        flatCols.push(c);
        flatHeaders.push(g.headers[i]);
        // <condition>_<replicate>; the condition may hold underscores (Ct_a_1)
        var mm = (g.headers[i] || '').match(/^(.+)_[0-9]+$/);
        if (mm) colCondMap[c] = mm[1];
      });
    });
    flatCols.forEach(function (c, i) { ws.getColumn(i + 1).width = x.widths[i]; });

    // Row 1: group headers (merged)
    var colOffset = 1;
    x.groups.forEach(function (g) {
      if (!g.columns || g.columns.length === 0) return;
      var startCol = colOffset;
      var endCol = colOffset + g.columns.length - 1;
      var cell = ws.getCell(1, startCol);
      cell.value = g.name;
      if (endCol > startCol) ws.mergeCells(1, startCol, 1, endCol);
      cell.font = { bold: true, color: { argb: 'FFFFFFFF' }, size: 12 };
      cell.fill = solid(x.headerFill);
      cell.alignment = { vertical: 'middle', horizontal: 'center' };
      colOffset = endCol + 1;
    });
    ws.getRow(1).height = 28;

    // Row 2: per-column headers, tinted by condition for the sample columns
    flatCols.forEach(function (c, i) {
      var cell = ws.getCell(2, i + 1);
      cell.value = flatHeaders[i];
      var cond = colCondMap[c];
      cell.fill = solid(cond && x.palette[cond] ? hex2argb(x.palette[cond]) : 'FF595959');
      cell.font = { bold: true, color: { argb: 'FFFFFFFF' }, size: 11 };
      cell.alignment = { vertical: 'middle', horizontal: 'center' };
    });
    ws.getRow(2).height = 24;

    // Data rows from row 3
    var nextRow = 3;
    rows.forEach(function (row) {
      var excelRow = ws.getRow(nextRow);
      flatCols.forEach(function (c, i) {
        var v = row[c];
        if (v != null && v !== '') {
          var num = text.indexOf(c) === -1 ? toNumber(v) : null;
          excelRow.getCell(i + 1).value = num !== null ? num : v;
        }
      });
      excelRow.commit();
      nextRow++;
    });

    for (var r = 1; r < nextRow; r++) {
      var used = ws.getRow(r);
      for (var c = 1; c <= flatCols.length; c++) used.getCell(c).border = thinBorder();
    }
    // Freeze both header rows and the first four columns
    ws.views = [{ state: 'frozen', xSplit: Math.min(4, flatCols.length), ySplit: 2 }];
  }

  // One row per sample, tinted by condition, with a strong colour on the
  // Condition and Coding cells (sample metadata)
  function sheetTinted(ws, rows, x) {
    var flatCols = x.columns;
    flatCols.forEach(function (c, i) { ws.getColumn(i + 1).width = x.widths[i]; });

    flatCols.forEach(function (c, i) {
      var cell = ws.getCell(1, i + 1);
      cell.value = c;
      cell.fill = solid(x.headerFill);
      cell.font = { bold: true, color: { argb: 'FFFFFFFF' }, size: 11 };
      cell.alignment = { vertical: 'middle', horizontal: 'center' };
    });
    ws.getRow(1).height = 28;

    var nextRow = 2;
    rows.forEach(function (row) {
      var cond = String(row['Condition'] || '');
      var tintArgb = x.tint[cond] ? hex2argb(x.tint[cond]) : 'FFFFFFFF';
      var strongArgb = x.palette[cond] ? hex2argb(x.palette[cond]) : 'FF808080';
      var excelRow = ws.getRow(nextRow);
      flatCols.forEach(function (c, i) {
        var cell = excelRow.getCell(i + 1);
        var numeric = x.numeric.indexOf(c) !== -1;
        var v = row[c];
        if (v != null && v !== '') {
          if (numeric) {
            var num = parseFloat(v);
            cell.value = !isNaN(num) ? num : v;
            cell.numFmt = '#,##0';
          } else {
            cell.value = v;
          }
        }
        if (x.strong.indexOf(c) !== -1) {
          cell.fill = solid(strongArgb);
          cell.font = { bold: true, color: { argb: 'FFFFFFFF' } };
          cell.alignment = { vertical: 'middle', horizontal: 'center' };
        } else {
          cell.fill = solid(tintArgb);
          cell.alignment = { vertical: 'middle',
                             horizontal: i === 0 ? 'left' : (numeric ? 'right' : 'center') };
        }
      });
      excelRow.commit();
      nextRow++;
    });

    for (var r = 1; r < nextRow; r++) {
      var used = ws.getRow(r);
      for (var c = 1; c <= flatCols.length; c++) used.getCell(c).border = thinBorder('FFBFBFBF');
    }
    ws.views = [{ state: 'frozen', xSplit: 1, ySplit: 1 }];
  }

  // --- Public interface -----------------------------------------------------

  window.NadiaTables = {
    register: function (id, settings) {
      settings.hidden = {};
      settings.filtersVisible = false;
      tables[id] = settings;
    },

    toggleFilters: function (id) {
      var t = tables[id];
      var container = root(id).querySelector('.rl-filters-container');
      var btn = root(id).querySelector('.rl-btn-toggle-filters');
      if (!container) return;
      var s = container.style;
      if (t.filtersVisible) {
        s.maxHeight = '0'; s.opacity = '0'; s.marginBottom = '0'; s.padding = '0';
        s.borderWidth = '0'; s.boxShadow = 'none'; s.overflow = 'hidden';
        if (btn) btn.classList.add('filters-hidden');
      } else {
        s.maxHeight = '500px'; s.opacity = '1'; s.marginBottom = '1rem'; s.padding = '1.5rem';
        s.borderWidth = '1px'; s.boxShadow = '0 2px 4px rgba(0, 0, 0, 0.05)'; s.overflow = 'visible';
        if (btn) btn.classList.remove('filters-hidden');
      }
      t.filtersVisible = !t.filtersVisible;
    },

    setFilter: function (id, column, value) {
      Reactable.setFilter(id, column, value || undefined);
    },

    toggleCondition: function (id, cond) {
      var t = tables[id];
      t.hidden[cond] = !t.hidden[cond];
      chips(id).forEach(function (chip) {
        if (chip.getAttribute('data-cond') === cond) chip.classList.toggle('off', t.hidden[cond]);
      });
      applyConditions(id);
    },

    clearFilters: function (id) {
      var t = tables[id];
      var r = root(id);
      r.querySelectorAll('.rl-filter-item .selectized').forEach(function (sel) {
        if (sel.selectize) sel.selectize.clear();
      });
      r.querySelectorAll('.rl-filters-container .rl-numeric-input').forEach(function (inp) {
        inp.value = '';
        var col = inp.getAttribute('data-column');
        if (col) Reactable.setFilter(id, col, undefined);
      });
      if (t.conditionMode) {
        Object.keys(t.hidden).forEach(function (c) { t.hidden[c] = false; });
        chips(id).forEach(function (chip) { chip.classList.remove('off'); });
        applyConditions(id);
      }
      var search = r.querySelector('.rl-search-input');
      if (search) {
        search.value = '';
        Reactable.setSearch(id, '');
      }
    },

    exportExcel: async function (id) {
      var x = tables[id].excel;
      try {
        var tsv = Reactable.getDataCSV(id, { sep: '\t' });
        var parsed = Papa.parse(tsv, { header: true, delimiter: '\t', skipEmptyLines: true });
        var wb = new ExcelJS.Workbook();
        var ws = wb.addWorksheet(x.sheet);
        if (x.layout === 'flat') sheetFlat(ws, parsed.data, parsed.meta.fields, x);
        else if (x.layout === 'grouped') sheetGrouped(ws, parsed.data, x);
        else sheetTinted(ws, parsed.data, x);

        var buffer = await wb.xlsx.writeBuffer();
        var blob = new Blob([buffer], { type: 'application/vnd.openxmlformats-officedocument.spreadsheetml.sheet' });
        var url = window.URL.createObjectURL(blob);
        var a = document.createElement('a');
        a.href = url;
        a.download = x.file + '_' + new Date().toISOString().split('T')[0] + '.xlsx';
        a.click();
        window.URL.revokeObjectURL(url);
      } catch (e) {
        console.error('Error while exporting:', e);
      }
    }
  };
})();
