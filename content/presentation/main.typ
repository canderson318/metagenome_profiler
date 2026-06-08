
#import "lib.typ": *

#let theme_col= rgb("5a0d87")
#let foot_col= red

#show footnote.entry: set text(foot_col,size:.8em)

#show footnote: it => (
    context {
      text(fill: foot_col, size: 1.2em)[#it]
    }
  )

#show figure.caption: emph

// Project configuration
#show: typslides.with(
  ratio: "16-9",
  theme: theme_col,
  font: "Fira Sans",
  font-size: 20pt,
  link-style: "color",
  show-progress: true,
)

#front-slide(
  title: "Metagenome Profiler",
  authors: "Christian Anderson",
  info: today.display( "[day] [month repr:long] [year]"),
)

#table-of-contents()


#title-slide[Background]



#slide(title: "Problem", outlined: false)[
  Hypothetical
]

#slide(outlined: false)[
  #framed(title: "Problem")[PROBLEM]
]

#title-slide[Approach]

#slide(title: "Algorithm outline", outlined: false)[
]

#title-slide[Results]

#slide(title: "Results", outlined: false)[
  #image("../../out/figs/")
]


#slide(title: "Findings", outlined: false)[
]


#slide(title: "Limitations & future directions", outlined: false)[
]


// ─────────────────────────────────────────────
// BIBLIOGRAPHY
// ─────────────────────────────────────────────
#let bib = bibliography("bibliography.bib")
#bibliography-slide(bib, extra: align(bottom)[ https://github.com/canderson318/research-methodology-programming-assignment ])
