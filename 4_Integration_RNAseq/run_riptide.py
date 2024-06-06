import riptide
import cobra
import pandas as pd

#model = cobra.io.read_sbml_model('/home/mferreira/Mauricio/Biologia/Doutorado/3.Kmax_ethanol/2_ecGEM_reconstruction/models/Kluyveromyces_marxianus-GEM.xml')
model = cobra.io.load_matlab_model('ecKmarx.mat')

for reaction in model.exchanges:
    reaction.lower_bound = -1000
    reaction.upper_bound = 1000

model.objective = 'r_1913'

Diniz2017_0h = riptide.read_transcription_file('Diniz2017_normcounts_0h.csv')
Diniz2017_1h = riptide.read_transcription_file('Diniz2017_normcounts_1h.csv')
Diniz2017_4h = riptide.read_transcription_file('Diniz2017_normcounts_4h.csv')

Mo2019_0KM = riptide.read_transcription_file('Mo2019_normcounts_0KM.csv')
Mo2019_4KM = riptide.read_transcription_file('Mo2019_normcounts_4KM.csv')
Mo2019_6KM = riptide.read_transcription_file('Mo2019_normcounts_6KM.csv')

Mo2019_0100d = riptide.read_transcription_file('Mo2019_normcounts_0100d.csv')
Mo2019_4100d = riptide.read_transcription_file('Mo2019_normcounts_4100d.csv')
Mo2019_6100d = riptide.read_transcription_file('Mo2019_normcounts_6100d.csv')

riptide_Diniz2017_0h  = riptide.contextualize(model=model, transcriptome=Diniz2017_0h)
riptide_Diniz2017_1h  = riptide.contextualize(model=model, transcriptome=Diniz2017_1h)
riptide_Diniz2017_4h  = riptide.contextualize(model=model, transcriptome=Diniz2017_4h)

riptide_Mo2019_0KM = riptide.contextualize(model=model, transcriptome=Mo2019_0KM)
riptide_Mo2019_4KM = riptide.contextualize(model=model, transcriptome=Mo2019_4KM)
riptide_Mo2019_6KM = riptide.contextualize(model=model, transcriptome=Mo2019_6KM)

riptide_Mo2019_0100d = riptide.contextualize(model=model, transcriptome=Mo2019_0100d)
riptide_Mo2019_4100d = riptide.contextualize(model=model, transcriptome=Mo2019_4100d)
riptide_Mo2019_6100d = riptide.contextualize(model=model, transcriptome=Mo2019_6100d)

riptide.save_output(riptide_obj=riptide_Diniz2017_0h, path='Diniz2017_0h', file_type="sbml")
riptide.save_output(riptide_obj=riptide_Diniz2017_1h, path='Diniz2017_1h', file_type="sbml")
riptide.save_output(riptide_obj=riptide_Diniz2017_4h, path='Diniz2017_4h', file_type="sbml")

riptide.save_output(riptide_obj=riptide_Mo2019_0KM, path='Mo2019_0KM', file_type="sbml")
riptide.save_output(riptide_obj=riptide_Mo2019_4KM, path='Mo2019_4KM', file_type="sbml")
riptide.save_output(riptide_obj=riptide_Mo2019_6KM, path='Mo2019_6KM', file_type="sbml")

riptide.save_output(riptide_obj=riptide_Mo2019_0100d, path='Mo2019_0100d', file_type="sbml")
riptide.save_output(riptide_obj=riptide_Mo2019_4100d, path='Mo2019_4100d', file_type="sbml")
riptide.save_output(riptide_obj=riptide_Mo2019_6100d, path='Mo2019_6100d', file_type="sbml")
