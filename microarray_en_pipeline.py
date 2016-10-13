# -*- coding: utf-8 -*-

import sys
from collections import defaultdict
from plate_well_info import getLevel2info
import numpy as np
from multiprocessing import Pool
from collections import Counter
import random
import pandas as pd

reload(sys)
sys.setdefaultencoding('UTF8')
sys.path.insert(0, '/Users/daniel/repos/python-client')  # for paginated get_file_children

from genestack_client import (FilesUtil, Application, GenestackException,
                              get_connection, SpecialFolders,
                              CLApplication, Metainfo)


class L1000NormalisationApplication(CLApplication):
    APPLICATION_ID = "genestack/l1000-normalisation"


HG_U133_ANNOTATION = "GSF12957143"  # int-dev
L1000_ANNOTATION = "GSF3659936"  # dotorg


class ExpressionNavigatorForMicroarraysApplication(Application):
    APPLICATION_ID = "genestack/expressionNavigator-microarrays"

    def create_file(self, parameters):
        return self.invoke('createFile', parameters)


print "Connecting to Genestack..."
connection = get_connection()
fu = FilesUtil(connection)
en = ExpressionNavigatorForMicroarraysApplication(connection)
norm_app = L1000NormalisationApplication(connection)
created_files = fu.get_special_folder(SpecialFolders.CREATED)

RAW_MICROARRAY_FOLDERS = ['GSF1943722']
                          #'GSF1952703', 'GSF1913719', 'GSF1933721', 'GSF1923720', 'GSF1893717', 'GSF1859015',
                          #'GSF1903718', 'GSF1943722', 'GSF1852785', 'GSF1883716', 'GSF1869016', 'GSF1962706']


PERTURBATION_TYPE = "sourceData:geo.Sample_treatment_protocol_ch1.SM_Pert_Type"
COMPOUND_PERTURBATION = "trt_cp"
TIME = "sourceData:geo.Sample_treatment_protocol_ch1.SM_Time"
TIME_UNIT = "sourceData:geo.Sample_treatment_protocol_ch1.SM_Time_Unit"
DOSE = "sourceData:geo.Sample_treatment_protocol_ch1.SM_Dose"
DOSE_UNIT = "sourceData:geo.Sample_treatment_protocol_ch1.SM_Dose_Unit"
COMPOUND = "sourceData:geo.Sample_treatment_protocol_ch1.SM_Name"
CELL_LINE = "sourceData:geo.Sample_characteristics_ch1.cl_name"

KEYS = [Metainfo.NAME, PERTURBATION_TYPE, TIME, CELL_LINE, TIME_UNIT, DOSE, DOSE_UNIT, COMPOUND]

plateInfo = getLevel2info()
n = len(plateInfo.id)
plateDict = {}
for x in range(n):
    plateDict[plateInfo.id[x]] = plateInfo.det_plate[x]

def flatten(l):
    return flatten(l[0]) + (flatten(l[1:]) if len(l) > 1 else []) if type(l) is list else [l]


dmsoAdded = {}

def group_DMSOs_for_analysis(infos):
    print "Collecting DMSOs.."
    DMSOs = defaultdict(set)

    for accession, info in infos.iteritems():
        info = infos[accession]
        compoundName = info[COMPOUND]
        plate = plateDict[info['genestack:name']].split('_')[0]
        if compoundName.lower() == "dmso":
            simpleSig = "{0} - {1} - {2}{3}".format(info[CELL_LINE],
                                              plate,
                                              info[TIME], info[TIME_UNIT])
            DMSOs[simpleSig].add(accession)
        signature = "{0} - {1} - {2} - {3}{4}".format(info[CELL_LINE], compoundName, plate, info[TIME], info[TIME_UNIT])
        dmsoAdded[signature] = False
    return DMSOs


print "Collecting microarray files..."
#all_marray_infos = {}
all_marray_infos = defaultdict(dict)

dfMeta = pd.read_csv('lincs_metadata.csv', sep='\t')
for x in range(len(dfMeta)):
    if x % 5000 == 0:
        print x, "out of", n
    all_marray_infos[dfMeta.iloc[x, 0]] = dict(dfMeta.iloc[x][KEYS])

    '''
    currentid = dfMeta.iloc[x,0]
    all_marray_infos.setdefault(currentid, [])
    all_marray_infos[currentid] = {}
    for k in KEYS:
        all_marray_infos[currentid][k] = dfMeta.ix[x][k]
    '''

'''

for fd in RAW_MICROARRAY_FOLDERS:
    marray_files = fu.get_file_children(fd)#, show_large_folder_progress=True)
    infos = fu.get_metainfo_values_as_strings(marray_files, KEYS)
    all_marray_infos.update(infos)
'''
# with open('lincs_metadata.csv', 'w') as out:
#     out.write('\t'.join(['accession'] + KEYS) + '\n')
#     for acc, dic in all_marray_infos.iteritems():
#         out.write('\t'.join([acc] + [dic.get(k) if dic.get(k) is not None else '-' for k in KEYS]) + inf2[acc]+'\n')

# sys.exit(0)
DMSOs = group_DMSOs_for_analysis(all_marray_infos)
print "Collected DMSO data: "
print [[k, len(DMSOs[k])] for k in DMSOs.keys()]

def group_samples_for_analysis(infos):
    groups = defaultdict(set)

    print "Collecting compounds..."

    for accession, info in infos.iteritems():
        info = infos[accession]
        compoundName = info[COMPOUND]
        plate = plateDict[info['genestack:name']].split('_')[0]
        if compoundName.lower() != "dmso":
            signature = "{0} - {1} - {2} - {3}{4}".format(info[CELL_LINE], compoundName, plate, info[TIME], info[TIME_UNIT])
            simpleSig = "{0} - {1} - {2}{3}".format(info[CELL_LINE], plate, info[TIME], info[TIME_UNIT])
            #print signature, len(DMSOs[simpleSig]), dmsoAdded[signature]
            if simpleSig in DMSOs:
                #if dmsoAdded == False:
                #    if len(DMSOs[simpleSig])>3:
                #        groups[signature].update(random.sample(DMSOs[simpleSig], 3))
                #        dmsoAdded[signature] = True
                #    else:
                groups[signature].update(DMSOs[simpleSig])
                #        dmsoAdded[signature] = True
            groups[signature].add(accession)
    return groups

grouped = group_samples_for_analysis(all_marray_infos)
print "Compounds groups have been successfully collected!"
print [[k, len(grouped[k])] for k in grouped.keys()[0:25]]

print "\nCollecting contrast info from source files..."

contrast_values = {ma_file: all_marray_infos[ma_file].get(DOSE) + all_marray_infos[ma_file][DOSE_UNIT]
                   if all_marray_infos[ma_file][DOSE]!='null' else '0'
                   for ma_file in all_marray_infos}
### DOSE edited for GSF1968849, GSF1967775, GSF1968848
links = {}
norm_folder = fu.create_folder("Normalised files for L1000 data")
norm_files_to_sources = {}
for i, (name, files) in enumerate(grouped.iteritems(), 1):
    files = list(files)
    sys.stdout.write("\rCreating normalised microarray files (%d/%d)..." % (i, len(grouped)))
    sys.stdout.flush()
    if len(files) == 1:
        continue
    norm_file = norm_app.create_file(files, name=name)
    norm_files_to_sources[norm_file] = files
    links[norm_file] = [norm_folder]
    if i > 100:
        break
fu.link_files(links)
fu.unlink_files({norm_file: [created_files] for norm_file in links})

# we need a contrast_values dict of microarray file accession -> dose ("25 mg")

created_files = fu.get_special_folder(SpecialFolders.CREATED)
parent = fu.create_folder("Differential Expression For L1000", created_files)
file_links = {}
countNoDose = 0

for i, norm_file in enumerate(norm_files_to_sources,1):
    try:
        sys.stdout.write("\rCreating EN files ({1}/{2})...".format(norm_file, i, len(norm_files_to_sources)))
        sys.stdout.flush()

        source_files = norm_files_to_sources[norm_file]
        unique_contrasts = list(set(contrast_values.get(source) for source in source_files))
        # group_mappings is a dict { unique_contrast_value -> index }
        group_mappings = {value: i for i, value in enumerate(unique_contrasts, 1)}
        group_names = ["Ungrouped"] + ["Group {0}".format(contrast) for contrast in unique_contrasts]
        program_options = {'microarrayAnnotationSource': {
            'type': 'source',
            'value': L1000_ANNOTATION
        }}

        # TODO find control groups for L1000
        # if there is a dose 'null', use it as control group
        control_index = None
        for k, contrast in enumerate(unique_contrasts, 1):
            if contrast == '0':
                control_index = k
        if control_index is not None:
            groupCounts = Counter([group_mappings[contrast_values[source]] for source in source_files])
            if all(v>1 for v in groupCounts.values()):
                program_options['controlGroupOption'] = {'value': str(control_index)}
            else:
                countNoDose += 1
                fu.unlink_file(norm_file,  norm_folder)
                continue
        params = {
            'organism': "Homo sapiens",
            'normalisedInputMode': True,
            'accessionList': source_files,
            'groupIdList': map(str, [group_mappings[contrast_values[source]] for source in source_files]),
            'groupsNameList': group_names,
            'groupsDescriptionList': group_names,
            'programOptions': program_options,
            'sourcesAccessions': [norm_file]
        }
        en_file = en.create_file(params)
        file_links[en_file] = [parent]
    except GenestackException as e:
        print "Received Genestack exception: {0}".format(e.message)
print "There were {0} out of {1} files without NULL dose! ".format(countNoDose, len(norm_files_to_sources))
print "\nMoving files..."
fu.link_files(file_links)
fu.unlink_files({key: [created_files] for key in file_links})