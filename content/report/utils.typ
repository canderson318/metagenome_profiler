
////////////////////////
// Function Definitions
////////////////////////

//************************************************************************\\
// Creates a filled rect with text
// Make a Todo note
//
// Arguments: 
//   body: Text to show
//   color: bg fill, (default: yellow)
//   title: rect label (default: TODO)
//
// Returns: a colored box with body text labeled TODO
#let TODO(body, color: yellow, title: "TODO") = {
  rect(
    width: 100%,
    radius: 3pt,
    stroke: 0.5pt,
    fill: color,
  )[
    #text(weight: 700)[#title]: #body
  ]
}

//************************************************************************\\

#let today = datetime.today()

//************************************************************************\\
