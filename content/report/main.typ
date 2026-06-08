
#import "utils.typ": *
#import "@preview/zebraw:0.6.1": *
#show: zebraw

////////////////////
// Document Settings
////////////////////
#set page(margin: .5in, width: 8.5in, height: 11in, numbering: "1")
#let fntsz = 11pt
#let after = fntsz * 1.5
#set text(font: "Georgia", size: fntsz)
#set par(leading: fntsz, spacing: after)
#set heading(numbering: "1.")

////////////////
// Set Defaults
////////////////

#show math.equation: set text(size: 12pt)

#show title: it => {
  set align(center + horizon)
  set text(size: 20pt, weight: "bold")
  v(-5em)
  it
}

#show outline: it => {
  show heading: set text(size: 15pt)
  show outline.entry: it => pad(left: 1em, right: 0em, it)
  it
}

#show link: it =>{
    set text(blue)
    it
}

#show ref: it =>{
    set text(style:"italic")
    it
}

#show heading: it=>{
    block(below:after)[#it]
}

#show figure.where(kind: "code"): it => {
    // let col = rgb("#2746e397")
    // block(fill: none, stroke: col, inset: 8pt,below:5pt, radius: 20pt, width: 100%)[
    //     #it.body
    // ]
    block[
        #it.body
    ]
    v(-1em)
    it.caption
}

#set heading(numbering:none)

#let theme_col= rgb("#8c2fc185")

#show "XXX": text.with(fill:red)

#show "GH": box[
    #box(
        link("https://github.com/canderson318/research-methodology-programming-assignment", image("../github.svg", height: .8em))
    )
]

//////////////////////
//////////////////////
//////// BODY ////////
//////////////////////
//////////////////////

// Written report (25%): This maximum 4-page (11pt arial font, ½ inch margins) technical document with an introduction, methods, and results sections. The introduction gives any necessary background and describes and justifies your strategy. The method section defines the input, presents the algorithm, and details the expected output. The results contain experiments, results, and discussion. You should cite the appropriate scientific literature where appropriate. The list of citations does not count against the page limit.

#title[Metatranscriptomic Virome Characterization of _Rhinolophus malayanus_ Bat Feces]
#align(center)[By: Christian Anderson\ #today.display("[day] [month repr:long] [year]")]


#align(bottom, outline())
#pagebreak()


= Introduction
// - assigned problem
// - Overview 
// - Computational translation
// - Biological Relevance
// - importance 
// - example 


= Methods

= Results

= Discussion   


= AI use statement

