# BayesPower 1.0.4

* Updated documentation for all exported functions.


# BayesPower 1.0.3

* Simplified S3 methods:
  * `plot.BFpower`
  * `print.BFpower`
  * `print.BFvalue`
* Renamed arguments in `BFpower.*` and `BF10*` functions:
  * `hypothesis` → `alternative`
  * `e` → `ROPE`
  * `model` → `prior_analysis`
  * `model_d` → `prior_design`
  * `D` → `threshold`
  * `direct` → `type_rate`
* Simplified `BFpower.*` and `BF10*` functions by reducing the number of required arguments.
* **Breaking change:** Argument and function renaming affect backward compatibility.


# BayesPower 1.0.2

* Added S3 `plot` and `print` methods for Bayes factor and power objects.
* Expanded and improved documentation for all exported functions.
* Added input validation checks for `BFpower.*` and `BF10*` functions.
* Renamed functions for consistency:
  * `BF10.t.test.one_sample` → `BF10.ttest.OneSample`
  * `BF10.t.test.two_sample` → `BF10.ttest.TwoSample`
  * `BFpower.t.test_one_sample` → `BFpower.ttest.OneSample`
  * `BFpower.t.test_two_sample` → `BFpower.ttest.TwoSample`
  * `BFpower.f` → `BFpower.f.test`
* **Breaking change:** Function renaming affect backward compatibility.


# BayesPower 1.0.1

* Initial release.
* Added Shiny-related function:
  * `BayesPower_BayesFactor`
* Implemented Bayes factor functions:
  * `BF10.bin.test`
  * `BF10.cor`
  * `BF10.f.test`
  * `BF10.props`
  * `BF10.t.test.one_sample`
  * `BF10.t.test.two_sample`
* Implemented power analysis functions:
  * `BFpower.bin`
  * `BFpower.cor`
  * `BFpower.f`
  * `BFpower.props`
  * `BFpower.t.test_one_sample`
  * `BFpower.t.test_two_sample`
