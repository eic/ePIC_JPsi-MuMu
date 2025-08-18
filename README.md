# ePIC_JPsi
Repository containing some simple analysis scripts as well as instructions for producing the necessary files required for analysis.

## lAger Instructions
Install version [3.61](https://eicweb.phy.anl.gov/monte_carlo/lager/-/tree/3.6.1?ref_type=tags)

Follow the install instructions in the README of that repo.

### 10x130
Run with the [.json file](https://github.com/smithalex0024/ePIC_JPsi/blob/main/lAger/10x130ep_00mrad/jpsi-10on130.disp-jpsi-00-muon.json) included in this repository:
`lager -c jpsi-10on130.disp-jpsi-00-muon.json -r 1 -o <OUTPUT_DIRECTORY>`

Process the lAger output through the afterburner using the config **ip6_hiacc_100x10**

```
abconv -p ip6_hiacc_100x10 <OUTPUT_DIRECTORY>/lager-vmp-00mrad.jpsi-10on130.4pi.disp-jpsi-00-muon.run00001-lumi1.hepmc -o <OUTPUT_DIRECTORY>/ab_output-00mrad.jpsi-10on130.4pi.disp-jpsi-00-muon.run00001-lumi10
```

### 10x250
Run with the [.json file](https://github.com/smithalex0024/ePIC_JPsi/blob/main/lAger/10x250ep_00mrad/jpsi-10on250.disp-jpsi-00-muon.json) included in this repository:
`lager -c jpsi-10on250.disp-jpsi-00-muon.json -r 1 -o <OUTPUT_DIRECTORY>`

Process the lAger output through the afterburner using the config **ip6_hiacc_275x10**

```
abconv -p ip6_hiacc_275x10 <OUTPUT_DIRECTORY>/lager-vmp-00mrad.jpsi-10on250.4pi.disp-jpsi-00-muon.run00001-lumi1.hepmc -o <OUTPUT_DIRECTORY>/ab_output-00mrad.jpsi-10on250.4pi.disp-jpsi-00-muon.run00001-lumi10
```
