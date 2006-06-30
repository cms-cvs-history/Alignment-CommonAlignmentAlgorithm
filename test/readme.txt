

- tenmuons.gen_sim_digi.cfg

Generate, simulate and digitize particle gun events with
10 muons per event, flat in eta [-2.5,2.5] and pt [5,50] GeV

- trackreco.cfg

Reads digi file, performs CKF track reconstruction,
and writes new reco file

- minireco.cfg

Reads reco file, writes minireco file (filtering tracks with pt>10 GeV
and putting them into TkAlDST)

- trackrefit.cfg

Reads minireco file and refits the tracks

- alignment.cfg

Reads minireco file and runs AlignmentProducer
