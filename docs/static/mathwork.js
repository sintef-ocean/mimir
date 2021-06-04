window.MathJax = {
  tex: {
    macros: {
      Z: "\\mathbb{Z}",
      dat: "\\dot",
      st: "\\operatorname{s.t}",
      col: "\\operatorname{col}",
      row: "\\operatorname{row}",
      vek: "\\operatorname{vec}",
      bdiag: "\\operatorname{bdiag}",
      diag: "\\operatorname{diag}",
      band: "\\operatorname{band}",
      repvert: "\\operatorname{repvert}",
      rephorz: "\\operatorname{rephorz}"
    },
    inlineMath: [["\\(", "\\)"]],
    displayMath: [["\\[", "\\]"]],
    processRefs: false,
    processEnvironments: false,
    autoload: {
      color: [],
      colorV2: ['color']
    },
    packages: {'[+]': ['noerrors']},
  },
  options: {
    ignoreHtmlClass: 'tex2jax_ignore',
    processHtmlClass: 'tex2jax_process'
  }
};
