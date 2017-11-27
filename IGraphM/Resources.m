
(* Icon for IGClusterData summary boxes *)
$igClusterDataIcon =
    With[
      {
        c1 = Hue[0.56, 0.4, 0.9],
        c2 = Hue[0.41, 0.4, 0.9],
        c3 = Hue[1, 0.4, 0.9]
      },
      Graphics[
        {
          {
            {
              EdgeForm @ Directive[AbsoluteThickness[4], c1],
              c1,
              FilledCurve[
                BSplineCurve[
                  {
                    {1.9643355098788902, 1.15632431724549},
                    {2.1129562552212175, 1.15632431724549},
                    {2.3449402850809546, 1.4263017932106137},
                    {2.3449402850809546, 1.5749225385529408},
                    {2.0781787052473715, 1.9227916943029346},
                    {1.929557959905044, 1.9227916943029346},
                    {1.680620534733193, 1.7808600917946806},
                    {1.680620534733193, 1.6322393464523535}
                  },
                  SplineClosed -> True
                ]
              ]
            },
            {
              EdgeForm @ Directive[AbsoluteThickness[4], c2],
              c2,
              FilledCurve[
                BSplineCurve[
                  {
                    {0.26155496896322245, 2.0650294545810315},
                    {-0.07431037267116364, 1.5683181368247658},
                    {-0.07431037267116364, 1.4196973914824387},
                    {0.07431037267116364, 1.4196973914824387},
                    {0.4986641259554429, 1.4884187429358733},
                    {0.4986641259554429, 1.6370394882782004},
                    {0.4101757143055497, 2.0650294545810315}
                  },
                  SplineClosed -> True
                ]
              ]
            },
            {
              EdgeForm @ Directive[AbsoluteThickness[4], c3],
              c3,
              FilledCurve[
                BSplineCurve[
                  {
                    {0.7949214963560174, 0.11573578412817837},
                    {1.3077575771793024, -0.07431037267116364},
                    {1.4563783225216296, -0.07431037267116364},
                    {1.4563783225216296, 0.07431037267116364},
                    {1.2144410769380514, 0.4453807951685744},
                    {1.0658203315957242, 0.4453807951685744},
                    {0.7949214963560174, 0.26435652947050564}
                  },
                  SplineClosed -> True
                ]
              ]
            }
          },
          {
            {
              GrayLevel[0.3], Opacity[0.5], AbsoluteThickness[1],
              Line @ {
                {2.0386458825500537, 1.2306346899166536},
                {2.270629912409791, 1.5006121658817773}
              },
              Line @ {
                {2.0386458825500537, 1.2306346899166536},
                {2.0038683325762077, 1.848481321631771}
              },
              Line @ {
                {2.0386458825500537, 1.2306346899166536},
                {1.7549309074043566, 1.706549719123517}
              },
              BezierCurve[
                {
                  {2.0386458825500537, 1.2306346899166536},
                  {1.8695996296395192, 1.0372301013701672},
                  {1.6232253667757388, 0.7538330208353062},
                  {1.3768511039119584, 0.4704359403004452},
                  {0.8692318690271811, 0.190046156799342}
                },
                SplineDegree -> 2
              ],
              BezierCurve[
                {
                  {2.0386458825500537, 1.2306346899166536},
                  {1.8695996296395192, 1.0372301013701672},
                  {1.6232253667757388, 0.7538330208353062},
                  {1.3768511039119584, 0.4704359403004452},
                  {1.382067949850466, 0.}
                },
                SplineDegree -> 2
              ],
              Line @ {
                {2.270629912409791, 1.5006121658817773},
                {2.0038683325762077, 1.848481321631771}
              },
              Line @ {
                {2.270629912409791, 1.5006121658817773},
                {1.7549309074043566, 1.706549719123517}
              },
              BezierCurve[
                {
                  {2.270629912409791, 1.5006121658817773},
                  {1.8695996296395192, 1.0372301013701672},
                  {1.6232253667757388, 0.7538330208353062},
                  {1.3768511039119584, 0.4704359403004452},
                  {1.1401307042668878, 0.37107042249741073}
                },
                SplineDegree -> 2
              ],
              Line @ {
                {2.0038683325762077, 1.848481321631771},
                {1.7549309074043566, 1.706549719123517}
              },
              BezierCurve[
                {
                  {2.0038683325762077, 1.848481321631771},
                  {1.5045771018575447, 1.7773176649728462},
                  {1.1297545837248073, 1.777119809568048},
                  {0.75493206559207, 1.7769219541632504},
                  {0.42435375328427927, 1.5627291156070369}
                },
                SplineDegree -> 2
              ],
              BezierCurve[
                {
                  {1.7549309074043566, 1.706549719123517},
                  {1.5045771018575447, 1.7773176649728462},
                  {1.1297545837248073, 1.777119809568048},
                  {0.75493206559207, 1.7769219541632504},
                  {0.3358653416343861, 1.990719081909868}
                },
                SplineDegree -> 2
              ],
              Line @ {
                {0.3358653416343861, 1.990719081909868},
                {0., 1.4940077641536023}
              },
              Line @ {
                {0.3358653416343861, 1.990719081909868},
                {0.42435375328427927, 1.5627291156070369}
              },
              Line @ {
                {0., 1.4940077641536023},
                {0.42435375328427927, 1.5627291156070369}
              },
              BezierCurve[
                {
                  {0., 1.4940077641536023},
                  {0.41030297914336333, 1.216415902322334},
                  {0.6084290816445871, 0.9044633647643481},
                  {0.8065551841458107, 0.5925108272063622},
                  {1.1401307042668878, 0.37107042249741073}
                },
                SplineDegree -> 2
              ],
              BezierCurve[
                {
                  {0.42435375328427927, 1.5627291156070369},
                  {0.41030297914336333, 1.216415902322334},
                  {0.6084290816445871, 0.9044633647643481},
                  {0.8065551841458107, 0.5925108272063622},
                  {0.8692318690271811, 0.190046156799342}
                },
                SplineDegree -> 2
              ],
              Line @ {
                {0.8692318690271811, 0.190046156799342},
                {1.382067949850466, 0.}
              },
              Line @ {
                {0.8692318690271811, 0.190046156799342},
                {1.1401307042668878, 0.37107042249741073}
              },
              Line @ {
                {1.382067949850466, 0.},
                {1.1401307042668878, 0.37107042249741073}
              }
            },
            {
              Black,
              Disk[{2.0386458825500537, 1.2306346899166536}, 0.1],
              Disk[{2.270629912409791, 1.5006121658817773}, 0.1],
              Disk[{2.0038683325762077, 1.848481321631771}, 0.1],
              Disk[{1.7549309074043566, 1.706549719123517}, 0.1],
              Disk[{0.3358653416343861, 1.990719081909868}, 0.1],
              Disk[{0., 1.4940077641536023}, 0.1],
              Disk[{0.42435375328427927, 1.5627291156070369}, 0.1],
              Disk[{0.8692318690271811, 0.190046156799342}, 0.1],
              Disk[{1.382067949850466, 0.}, 0.1],
              Disk[{1.1401307042668878, 0.37107042249741073}, 0.1]
            }
          }
        },
        Frame -> False,
        PlotRange -> {{-0.3, 2.6}, {-0.3, 2.3}},
        ImageSize -> Dynamic[{Automatic, 3.5 CurrentValue["FontCapHeight"] / AbsoluteCurrentValue[Magnification]}]
      ]
    ];
