###Code to study vertex reconstruction.

This assumes at least Root v6.22 and Cmake 3.15 but should work with slightly older versions, just make the change in the `CMakeLists.txt` file. To compile do the following,
```python
mkdir build
```
and set up using cmake
```python
cmake . -B build
```
then compile
```python
cd build && make
```
To run do
```python
./reco_vtx_study
```
with the output Root file `./build/out.root`

---

To add histograms open `hists.json` and follow the examples for either 1D
```json
    {
      "name" : "hEndXDiffEl",
      "type": "TH1D",
      "axes": "<Title>>;<Xaxis>>;<Yaxis>>",
      "bins": 50,
      "u_lim": 200,
      "l_lim": -300
    }
```
or 2D histograms
```json
    {
      "name" : "Example 2D Histo",
      "type": "TH2I",
      "axes": "<Title>>;<Xaxis>>;<Yaxis>>",
      "xbins": 10,
      "ybins": 10,
      "xu_lim": 10,
      "xl_lim": 0,
      "yu_lim": 10,
      "yl_lim": 0
    }
```
and add in the code as `Histogram_object.th1_hists["name_of_histogram"]->Fill( variable )`