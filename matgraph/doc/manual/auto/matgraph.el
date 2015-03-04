(TeX-add-style-hook "matgraph"
 (lambda ()
    (LaTeX-add-labels
     "sect:getting"
     "sect:using"
     "subsect:basic-principles"
     "subsect:getting-started"
     "subsection:declare"
     "subsect:inspectors"
     "subsect:graph-matrix"
     "subsect:invariants"
     "subsect:io"
     "subsect:large"
     "sect:documentation"
     "sect:perms"
     "sect:partition"
     "sect:under-the-hood"
     "inside-GM.g"
     "Q")
    (TeX-add-symbols
     "matlab"
     "ER"
     "RR"
     "ZZ"
     "oitem")
    (TeX-run-style-hooks
     "multicol"
     "graphicx"
     "superdate"
     "pslatex"
     "hyperref"
     "geometry"
     "latex2e"
     "amsart10"
     "amsart")))

