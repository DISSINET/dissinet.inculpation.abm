extensions [nw]

; In this model, there are two types of ties: connections (the input) and incriminations (the output)
undirected-link-breed [connections connection]
directed-link-breed [incriminations incrimination]

globals [
  ; CONNECTIONS
  count-conn         ; the number of connection ties (these ties are undirected) at the start
  count-conn-end     ; the number of connection ties at the end of the simulation
  clustering         ; clustering coefficient
  clustering-end     ; clustering coefficient at the end of the simulation
  components         ; number of components
  components-end     ; number of components at the end of the simulation
  degrees            ; degree distribution at the start
  degrees-end        ; degree distribution at the end of the simulation
  ; INCRIMINATIONS
  count-incr         ; the number of incrimination ties (these ties are directed)
  ties-incr          ; the incrimininations reported
  %deposed           ; percentage of agents deposed
  %incriminated      ; percentage of agents incriminated by at least another agent
  %incriminating     ; percentage of agents who incriminated at least another agent
  outdegrees-incr    ; outdegree distribution of incriminations
  indegrees-incr     ; indegree distribution of incriminations
]

turtles-own[
  deposed?           ; whether the agent has been deposed or not
  incriminated?      ; whether the agent has been incriminated by at least another agent
  incriminating?     ; whether the agent has incriminated at least another agent
  degree-conn        ; number of connection ties at the start
  degree-conn-end    ; number of connection ties at the end of the simulation
  outdegree-incr     ; number of incriminations sent
  indegree-incr      ; number of incriminations received
  betweenness        ; betweenness centrality
  temp_prob          ; auxiliary variable that is needed in the network generation process (if grow)
]

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; SETUP PROCEDURE
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
to setup
  clear-all
  reset-ticks

  ; two topologies availables: small-worlds and scale-free networks
  if topology = "small-world" [create-sw]
  if topology = "scale-free" [create-sf]
  if topology = "grow" [grow-network]

  ; These customise the agents and links at the start of the simulation
  customise-turtles
  customise-ties
  setup-globals
  update-plots
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; TEST PROCEDURE
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
to test
  ; the inquisitor can either pick deponents at random or snowball from previously incriminated individuals
  ifelse (snowballing? = true) [
    ; if there is no incriminated subject, they will pick one at random
    ifelse not any? turtles with [deposed? = false and incriminated? = true][
      ask one-of turtles with [deposed? = false] [depose]
    ][
      ask one-of turtles with [deposed? = false and incriminated? = true] [depose]
    ]
  ][
    ask one-of turtles with [deposed? = false] [depose]
  ]

  update-turtles
  update-ties
  update-globals

  print ties-incr ; print incriminations
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; GO PROCEDURE
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
to go
  ifelse any? turtles with [deposed? = false] [
    ifelse snowballing? [
      ifelse not any? turtles with [deposed? = false and incriminated? = true][
        ask one-of turtles with [deposed? = false] [depose]
      ][
        ask one-of turtles with [deposed? = false and incriminated? = true] [depose]
      ]
    ][
      ask one-of turtles with [deposed? = false] [depose]
    ]

    update-turtles
    update-ties
    update-globals
    tick
  ][
    ask turtles with [incriminated?][die]
    create-final-globals
    stop
  ]
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; NETWORK TOPOLOGIES
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
to create-sw
  nw:generate-watts-strogatz turtles connections N m rewiring-prob
  ; this is just for the layout on the screen
  layout-circle (sort turtles) max-pxcor - 5
  layout-circle (sort turtles with [who mod 2 = 0] ) max-pxcor - 10
  if rewiring-prob > 0 [
    repeat 5 [layout-spring turtles connections 0.2 5 1]
  ]
end

to create-sf
  nw:generate-preferential-attachment turtles connections N m
  ; this is just for the layout on the screen
  layout-circle (sort turtles) max-pxcor - 5
  layout-circle (sort turtles with [who mod 2 = 0] ) max-pxcor - 10
  repeat 5 [layout-spring turtles connections 0.2 5 1]
end

to grow-network
  create-turtles N
  ask turtles[
    setxy random-xcor random-ycor ; STEP 1: set agent's place on the environment at random
    let ego [who] of self
    let tie_count 0

    while [tie_count < m][
      let set_potential_interactors turtles with [(link-neighbor? turtle ego = false) and (who != [who] of turtle ego)]
      if count set_potential_interactors > 0[
        ask set_potential_interactors [set temp_prob exp( -1 * z * distance turtle ego)]
        let alter nobody
        let pick random-float (sum [temp_prob] of set_potential_interactors)
        let comp_low 0
        let comp_high 0
        ask set_potential_interactors
        [
          if alter = nobody
          [
            set comp_high comp_high + temp_prob
            ifelse pick >= comp_low and pick < comp_high
            [
              set alter who
            ]
            [
              set comp_low comp_high
            ]
          ]
        ]
        ask turtle alter [create-connection-with turtle ego]
        set tie_count tie_count + 1
      ]
    ]
  ]
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; MAIN PROCEDURES
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
to depose
  set deposed? true ; changed deposed? to true
  let i [who] of self ; ID of the deponent
  let potential-targets [who] of connection-neighbors ; list with the ID of connection-neighbours (it can be empty)

  if not empty? potential-targets[
    foreach potential-targets [ j ->
      ifelse (search-hubs? = true) [
        let x ([betweenness] of turtle j)
        let y random-normal x (x - x ^ 2) ; mean equal to between-centrality (x) and SD equal to (x - x ^ 2)
        if y < 0 [set y 0]
        if y > 1 [set y 1]
        if y > incrimination-threshold[
          ask turtle i [create-incrimination-to turtle j] ; create the incrimination
          set ties-incr lput (word i "-" j) ties-incr ; save the incrimination in the list in globals
          ask turtle j [set incriminated? true] ; change the stated of the target to incriminated
        ]
      ][
        let y random-float 1
        if y > incrimination-threshold[
          ask turtle i [create-incrimination-to turtle j] ; create the incrimination
          set ties-incr lput (word i "-" j) ties-incr ; save the incrimination in the list in globals
          ask turtle j [set incriminated? true] ; change the stated of the target to incriminated
        ]
      ]
    ]
  ]
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; CUSTOMISATION PROCEDURES
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
to customise-turtles
  ask turtles[
    set shape "circle"
    set size .6
    set color (grey + 2)
    set deposed? false
    set incriminated? false
    set incriminating? false
    set degree-conn count connection-neighbors
    set indegree-incr 0
    set outdegree-incr 0
    ; min-max feature scaling of betweenness centrality
    set betweenness ([nw:betweenness-centrality] of self - min [nw:betweenness-centrality] of turtles) /
    (max [nw:betweenness-centrality] of turtles - min [nw:betweenness-centrality] of turtles)
    ; if hide-labels? is true, agents will not show their ID number on the screen
    if (hide-labels? = false) [
      set label [who] of self
      set label-color yellow
    ]
  ]
end

to customise-ties
  ask connections [
    set color grey
    set thickness 0.1
    ; if hide-connections? is true, connections are not shown on the screen
    if (hide-connections? = true) [
      hide-link
    ]
  ]
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; UPDATE PROCEDURES
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
to setup-globals
  set count-conn count connections ; count the number of connections
  set count-conn-end "NA"
  set clustering mean [nw:clustering-coefficient] of turtles
  set clustering-end "NA"
  set components length nw:weak-component-clusters
  set components-end "NA"
  set degrees [degree-conn] of turtles
  set ties-incr [] ; incriminations is an empty list at the start
end

to update-turtles
  ask turtles[
    ; if deposed, agents change their shape to a triangle
    if deposed? = true[
      set shape "triangle"
      set size 1.5
    ]
    ; if incriminated, agents change their color to magenta
    if incriminated? = true[
      set color magenta
    ]

    set indegree-incr count in-incrimination-neighbors
    set outdegree-incr count out-incrimination-neighbors

    ; if somebody incriminated, set incriminating? to true
    if (outdegree-incr > 0) [
      set incriminating? true
    ]
  ]
end

to update-ties
  ask incriminations[
    set color (red + 1)
    set thickness 0.25
  ]
end

to update-globals
  set count-incr count incriminations ; count the number of incriminations
  set %deposed count turtles with [deposed?] / count turtles * 100
  set %incriminated count turtles with [incriminated?] / count turtles * 100
  set %incriminating count turtles with [incriminating?] / count turtles * 100
  set outdegrees-incr [outdegree-incr] of turtles
  set indegrees-incr [indegree-incr] of turtles
end

to create-final-globals
  set count-conn-end count connections
  if any? turtles [
    set clustering-end mean [nw:clustering-coefficient] of turtles
  ]
  set components-end length nw:weak-component-clusters
  ask turtles [
    set degree-conn-end count connection-neighbors
  ]
  set degrees-end [degree-conn-end] of turtles
  update-plots
end
@#$#@#$#@
GRAPHICS-WINDOW
189
10
602
424
-1
-1
5.0
1
11
1
1
1
0
0
0
1
-40
40
-40
40
1
1
1
ticks
30.0

CHOOSER
5
31
176
76
topology
topology
"small-world" "scale-free" "grow"
0

SLIDER
4
99
176
132
N
N
100
5000
300.0
100
1
NIL
HORIZONTAL

SLIDER
4
155
176
188
m
m
1
10
5.0
1
1
NIL
HORIZONTAL

SLIDER
4
193
176
226
rewiring-prob
rewiring-prob
0
1
0.23
0.01
1
NIL
HORIZONTAL

BUTTON
434
438
489
483
NIL
setup
NIL
1
T
OBSERVER
NIL
S
NIL
NIL
1

BUTTON
544
438
599
483
NIL
go
T
1
T
OBSERVER
NIL
G
NIL
NIL
1

TEXTBOX
10
81
160
99
Nodes
11
0.0
1

TEXTBOX
7
137
157
155
Connections created by node
11
0.0
1

BUTTON
489
438
544
483
NIL
test
NIL
1
T
OBSERVER
NIL
T
NIL
NIL
1

MONITOR
5
438
147
483
Connetions at the start
count-conn
17
1
11

MONITOR
609
556
706
601
Incriminations
count-incr
17
1
11

TEXTBOX
10
10
202
50
WORLD DIMENSIONS
15
105.0
1

TEXTBOX
8
304
227
322
BEHAVIOUR DIMENSIONS
15
105.0
1

SWITCH
5
230
176
263
hide-labels?
hide-labels?
0
1
-1000

SWITCH
5
267
176
300
hide-connections?
hide-connections?
1
1
-1000

PLOT
214
512
606
688
Degree distributions
degree
freq
0.0
20.0
0.0
10.0
true
true
"" ""
PENS
"Connections at start" 1.0 1 -7500403 true "" "histogram degrees"
"Incriminations sent" 1.0 1 -2674135 true "" "histogram outdegrees-incr"
"Incriminations received" 1.0 1 -13345367 true "" "histogram indegrees-incr"
"Connections at the end" 1.0 1 -16777216 true "" "histogram degrees-end"

MONITOR
609
512
706
557
% Deposed
%deposed
3
1
11

MONITOR
609
643
706
688
% Incriminated
%incriminated
3
1
11

MONITOR
609
600
706
645
% Incriminating
%incriminating
3
1
11

PLOT
8
512
208
688
% Incriminated
tick
percentage
0.0
10.0
0.0
100.0
true
false
"" ""
PENS
"default" 1.0 0 -2674135 true "" "plot count turtles with [incriminated?] / count turtles * 100"

MONITOR
147
438
282
483
Clustering at the start
clustering
3
1
11

MONITOR
8
692
150
737
Connections at the end
count-conn-end
0
1
11

TEXTBOX
12
489
162
508
OUTPUTS
15
105.0
1

MONITOR
150
692
279
737
Clustering at the end
clustering-end
3
1
11

MONITOR
282
438
421
483
Components at the start
components
0
1
11

MONITOR
279
692
422
737
Components at the end
components-end
0
1
11

SLIDER
4
397
177
430
incrimination-threshold
incrimination-threshold
0
1
0.92
0.01
1
NIL
HORIZONTAL

SWITCH
5
325
176
358
snowballing?
snowballing?
0
1
-1000

SWITCH
5
361
176
394
search-hubs?
search-hubs?
0
1
-1000

PLOT
422
592
606
737
Betweenness distribution
NIL
NIL
0.0
1.0
0.0
10.0
true
false
"" ""
PENS
"default" 0.01 1 -16777216 true "" "histogram [betweenness] of turtles"

INPUTBOX
605
10
655
70
z
10.0
1
0
Number

@#$#@#$#@
## WHAT IS IT?

(a general understanding of what the model is trying to show or explain)

## HOW IT WORKS

(what rules the agents use to create the overall behavior of the model)

## HOW TO USE IT

(how to use the model, including a description of each of the items in the Interface tab)

## THINGS TO NOTICE

(suggested things for the user to notice while running the model)

## THINGS TO TRY

(suggested things for the user to try to do (move sliders, switches, etc.) with the model)

## EXTENDING THE MODEL

(suggested things to add or change in the Code tab to make the model more complicated, detailed, accurate, etc.)

## NETLOGO FEATURES

(interesting or unusual features of NetLogo that the model uses, particularly in the Code tab; or where workarounds were needed for missing features)

## RELATED MODELS

(models in the NetLogo Models Library and elsewhere which are of related interest)

## CREDITS AND REFERENCES

(a reference to the model's URL on the web if it has one, as well as any other necessary credits, citations, and links)
@#$#@#$#@
default
true
0
Polygon -7500403 true true 150 5 40 250 150 205 260 250

airplane
true
0
Polygon -7500403 true true 150 0 135 15 120 60 120 105 15 165 15 195 120 180 135 240 105 270 120 285 150 270 180 285 210 270 165 240 180 180 285 195 285 165 180 105 180 60 165 15

arrow
true
0
Polygon -7500403 true true 150 0 0 150 105 150 105 293 195 293 195 150 300 150

box
false
0
Polygon -7500403 true true 150 285 285 225 285 75 150 135
Polygon -7500403 true true 150 135 15 75 150 15 285 75
Polygon -7500403 true true 15 75 15 225 150 285 150 135
Line -16777216 false 150 285 150 135
Line -16777216 false 150 135 15 75
Line -16777216 false 150 135 285 75

bug
true
0
Circle -7500403 true true 96 182 108
Circle -7500403 true true 110 127 80
Circle -7500403 true true 110 75 80
Line -7500403 true 150 100 80 30
Line -7500403 true 150 100 220 30

butterfly
true
0
Polygon -7500403 true true 150 165 209 199 225 225 225 255 195 270 165 255 150 240
Polygon -7500403 true true 150 165 89 198 75 225 75 255 105 270 135 255 150 240
Polygon -7500403 true true 139 148 100 105 55 90 25 90 10 105 10 135 25 180 40 195 85 194 139 163
Polygon -7500403 true true 162 150 200 105 245 90 275 90 290 105 290 135 275 180 260 195 215 195 162 165
Polygon -16777216 true false 150 255 135 225 120 150 135 120 150 105 165 120 180 150 165 225
Circle -16777216 true false 135 90 30
Line -16777216 false 150 105 195 60
Line -16777216 false 150 105 105 60

car
false
0
Polygon -7500403 true true 300 180 279 164 261 144 240 135 226 132 213 106 203 84 185 63 159 50 135 50 75 60 0 150 0 165 0 225 300 225 300 180
Circle -16777216 true false 180 180 90
Circle -16777216 true false 30 180 90
Polygon -16777216 true false 162 80 132 78 134 135 209 135 194 105 189 96 180 89
Circle -7500403 true true 47 195 58
Circle -7500403 true true 195 195 58

circle
false
0
Circle -7500403 true true 0 0 300

circle 2
false
0
Circle -7500403 true true 0 0 300
Circle -16777216 true false 30 30 240

cow
false
0
Polygon -7500403 true true 200 193 197 249 179 249 177 196 166 187 140 189 93 191 78 179 72 211 49 209 48 181 37 149 25 120 25 89 45 72 103 84 179 75 198 76 252 64 272 81 293 103 285 121 255 121 242 118 224 167
Polygon -7500403 true true 73 210 86 251 62 249 48 208
Polygon -7500403 true true 25 114 16 195 9 204 23 213 25 200 39 123

cylinder
false
0
Circle -7500403 true true 0 0 300

dot
false
0
Circle -7500403 true true 90 90 120

face happy
false
0
Circle -7500403 true true 8 8 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Polygon -16777216 true false 150 255 90 239 62 213 47 191 67 179 90 203 109 218 150 225 192 218 210 203 227 181 251 194 236 217 212 240

face neutral
false
0
Circle -7500403 true true 8 7 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Rectangle -16777216 true false 60 195 240 225

face sad
false
0
Circle -7500403 true true 8 8 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Polygon -16777216 true false 150 168 90 184 62 210 47 232 67 244 90 220 109 205 150 198 192 205 210 220 227 242 251 229 236 206 212 183

fish
false
0
Polygon -1 true false 44 131 21 87 15 86 0 120 15 150 0 180 13 214 20 212 45 166
Polygon -1 true false 135 195 119 235 95 218 76 210 46 204 60 165
Polygon -1 true false 75 45 83 77 71 103 86 114 166 78 135 60
Polygon -7500403 true true 30 136 151 77 226 81 280 119 292 146 292 160 287 170 270 195 195 210 151 212 30 166
Circle -16777216 true false 215 106 30

flag
false
0
Rectangle -7500403 true true 60 15 75 300
Polygon -7500403 true true 90 150 270 90 90 30
Line -7500403 true 75 135 90 135
Line -7500403 true 75 45 90 45

flower
false
0
Polygon -10899396 true false 135 120 165 165 180 210 180 240 150 300 165 300 195 240 195 195 165 135
Circle -7500403 true true 85 132 38
Circle -7500403 true true 130 147 38
Circle -7500403 true true 192 85 38
Circle -7500403 true true 85 40 38
Circle -7500403 true true 177 40 38
Circle -7500403 true true 177 132 38
Circle -7500403 true true 70 85 38
Circle -7500403 true true 130 25 38
Circle -7500403 true true 96 51 108
Circle -16777216 true false 113 68 74
Polygon -10899396 true false 189 233 219 188 249 173 279 188 234 218
Polygon -10899396 true false 180 255 150 210 105 210 75 240 135 240

house
false
0
Rectangle -7500403 true true 45 120 255 285
Rectangle -16777216 true false 120 210 180 285
Polygon -7500403 true true 15 120 150 15 285 120
Line -16777216 false 30 120 270 120

leaf
false
0
Polygon -7500403 true true 150 210 135 195 120 210 60 210 30 195 60 180 60 165 15 135 30 120 15 105 40 104 45 90 60 90 90 105 105 120 120 120 105 60 120 60 135 30 150 15 165 30 180 60 195 60 180 120 195 120 210 105 240 90 255 90 263 104 285 105 270 120 285 135 240 165 240 180 270 195 240 210 180 210 165 195
Polygon -7500403 true true 135 195 135 240 120 255 105 255 105 285 135 285 165 240 165 195

line
true
0
Line -7500403 true 150 0 150 300

line half
true
0
Line -7500403 true 150 0 150 150

pentagon
false
0
Polygon -7500403 true true 150 15 15 120 60 285 240 285 285 120

person
false
0
Circle -7500403 true true 110 5 80
Polygon -7500403 true true 105 90 120 195 90 285 105 300 135 300 150 225 165 300 195 300 210 285 180 195 195 90
Rectangle -7500403 true true 127 79 172 94
Polygon -7500403 true true 195 90 240 150 225 180 165 105
Polygon -7500403 true true 105 90 60 150 75 180 135 105

plant
false
0
Rectangle -7500403 true true 135 90 165 300
Polygon -7500403 true true 135 255 90 210 45 195 75 255 135 285
Polygon -7500403 true true 165 255 210 210 255 195 225 255 165 285
Polygon -7500403 true true 135 180 90 135 45 120 75 180 135 210
Polygon -7500403 true true 165 180 165 210 225 180 255 120 210 135
Polygon -7500403 true true 135 105 90 60 45 45 75 105 135 135
Polygon -7500403 true true 165 105 165 135 225 105 255 45 210 60
Polygon -7500403 true true 135 90 120 45 150 15 180 45 165 90

sheep
false
15
Circle -1 true true 203 65 88
Circle -1 true true 70 65 162
Circle -1 true true 150 105 120
Polygon -7500403 true false 218 120 240 165 255 165 278 120
Circle -7500403 true false 214 72 67
Rectangle -1 true true 164 223 179 298
Polygon -1 true true 45 285 30 285 30 240 15 195 45 210
Circle -1 true true 3 83 150
Rectangle -1 true true 65 221 80 296
Polygon -1 true true 195 285 210 285 210 240 240 210 195 210
Polygon -7500403 true false 276 85 285 105 302 99 294 83
Polygon -7500403 true false 219 85 210 105 193 99 201 83

square
false
0
Rectangle -7500403 true true 30 30 270 270

square 2
false
0
Rectangle -7500403 true true 30 30 270 270
Rectangle -16777216 true false 60 60 240 240

star
false
0
Polygon -7500403 true true 151 1 185 108 298 108 207 175 242 282 151 216 59 282 94 175 3 108 116 108

target
false
0
Circle -7500403 true true 0 0 300
Circle -16777216 true false 30 30 240
Circle -7500403 true true 60 60 180
Circle -16777216 true false 90 90 120
Circle -7500403 true true 120 120 60

tree
false
0
Circle -7500403 true true 118 3 94
Rectangle -6459832 true false 120 195 180 300
Circle -7500403 true true 65 21 108
Circle -7500403 true true 116 41 127
Circle -7500403 true true 45 90 120
Circle -7500403 true true 104 74 152

triangle
false
0
Polygon -7500403 true true 150 30 15 255 285 255

triangle 2
false
0
Polygon -7500403 true true 150 30 15 255 285 255
Polygon -16777216 true false 151 99 225 223 75 224

truck
false
0
Rectangle -7500403 true true 4 45 195 187
Polygon -7500403 true true 296 193 296 150 259 134 244 104 208 104 207 194
Rectangle -1 true false 195 60 195 105
Polygon -16777216 true false 238 112 252 141 219 141 218 112
Circle -16777216 true false 234 174 42
Rectangle -7500403 true true 181 185 214 194
Circle -16777216 true false 144 174 42
Circle -16777216 true false 24 174 42
Circle -7500403 false true 24 174 42
Circle -7500403 false true 144 174 42
Circle -7500403 false true 234 174 42

turtle
true
0
Polygon -10899396 true false 215 204 240 233 246 254 228 266 215 252 193 210
Polygon -10899396 true false 195 90 225 75 245 75 260 89 269 108 261 124 240 105 225 105 210 105
Polygon -10899396 true false 105 90 75 75 55 75 40 89 31 108 39 124 60 105 75 105 90 105
Polygon -10899396 true false 132 85 134 64 107 51 108 17 150 2 192 18 192 52 169 65 172 87
Polygon -10899396 true false 85 204 60 233 54 254 72 266 85 252 107 210
Polygon -7500403 true true 119 75 179 75 209 101 224 135 220 225 175 261 128 261 81 224 74 135 88 99

wheel
false
0
Circle -7500403 true true 3 3 294
Circle -16777216 true false 30 30 240
Line -7500403 true 150 285 150 15
Line -7500403 true 15 150 285 150
Circle -7500403 true true 120 120 60
Line -7500403 true 216 40 79 269
Line -7500403 true 40 84 269 221
Line -7500403 true 40 216 269 79
Line -7500403 true 84 40 221 269

wolf
false
0
Polygon -16777216 true false 253 133 245 131 245 133
Polygon -7500403 true true 2 194 13 197 30 191 38 193 38 205 20 226 20 257 27 265 38 266 40 260 31 253 31 230 60 206 68 198 75 209 66 228 65 243 82 261 84 268 100 267 103 261 77 239 79 231 100 207 98 196 119 201 143 202 160 195 166 210 172 213 173 238 167 251 160 248 154 265 169 264 178 247 186 240 198 260 200 271 217 271 219 262 207 258 195 230 192 198 210 184 227 164 242 144 259 145 284 151 277 141 293 140 299 134 297 127 273 119 270 105
Polygon -7500403 true true -1 195 14 180 36 166 40 153 53 140 82 131 134 133 159 126 188 115 227 108 236 102 238 98 268 86 269 92 281 87 269 103 269 113

x
false
0
Polygon -7500403 true true 270 75 225 30 30 225 75 270
Polygon -7500403 true true 30 75 75 30 270 225 225 270
@#$#@#$#@
NetLogo 6.2.2
@#$#@#$#@
@#$#@#$#@
@#$#@#$#@
<experiments>
  <experiment name="experiment" repetitions="50" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <metric>count-conn</metric>
    <metric>count-conn-end</metric>
    <metric>clustering</metric>
    <metric>clustering-end</metric>
    <metric>components</metric>
    <metric>components-end</metric>
    <metric>degrees</metric>
    <metric>degrees-end</metric>
    <metric>count-incr</metric>
    <metric>%deposed</metric>
    <metric>%incriminated</metric>
    <metric>%incriminating</metric>
    <metric>outdegrees-incr</metric>
    <metric>indegrees-incr</metric>
    <metric>ties-incr</metric>
    <enumeratedValueSet variable="topology">
      <value value="&quot;small-world&quot;"/>
      <value value="&quot;scale-free&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="N">
      <value value="150"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="m">
      <value value="1"/>
      <value value="2"/>
      <value value="3"/>
      <value value="4"/>
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="rewiring-prob">
      <value value="0.1"/>
      <value value="0.25"/>
      <value value="0.5"/>
      <value value="0.75"/>
      <value value="0.9"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="hide-labels?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="hide-connections?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="snowballing?">
      <value value="true"/>
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="search-hubs?">
      <value value="true"/>
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="incrimination-threshold">
      <value value="0.01"/>
      <value value="0.05"/>
      <value value="0.1"/>
      <value value="0.15"/>
      <value value="0.5"/>
      <value value="0.85"/>
      <value value="0.9"/>
      <value value="0.95"/>
      <value value="0.99"/>
    </enumeratedValueSet>
  </experiment>
</experiments>
@#$#@#$#@
@#$#@#$#@
default
0.0
-0.2 0 0.0 1.0
0.0 1 1.0 0.0
0.2 0 0.0 1.0
link direction
true
0
Line -7500403 true 150 150 120 180
Line -7500403 true 150 150 180 180
@#$#@#$#@
0
@#$#@#$#@
