import glob
import pandas as pd
import Bio.SeqUtils
from Bio import SeqIO


def extract_peaks_from_msd_file(file_list):
    '''extracts peaks from a list of .msd file generated with the software
    mmass (http://www.mmass.org/). For each input file the
    extracted information (m/z, intensity) is written to a tab-separated
    output file.

    Args:
        file_list (list): A list of .msd files
    '''
    for file in file_list:
        name = file.split('_')[-1].split('.')[0]
        to_write = list()
        with open(file) as in_file:
            for line in in_file:
                line = line.strip()
                if line.startswith('<peak mz="'):
                    line = line.split('"')
                    m_over_z = line[1]
                    intensity = line[3]
                    to_write.append((m_over_z, intensity))
        with open(name+'.txt', 'wt') as out_file:
            for item in to_write:
                out_file.write('\t'.join(item)+'\n')


class PeakCluster:
    def __init__(self, max_fraction, peak_list):
        self.max_fraction = max_fraction
        self.peak_list = peak_list
        self.name_list = list()


    def __lt__(self, other):
        self.max_int_value = self.peak_list[0].my_m_over_z
        other.max_int_value = other.peak_list[0].my_m_over_z
        return self.max_int_value < other.max_int_value


class Peak:
    def __init__(self, my_m_over_z, my_intensity, my_fraction):
        self.my_m_over_z = my_m_over_z
        self.my_intensity = my_intensity
        self.my_fraction = my_fraction


def compare_ms_spectra(peptide_file):
    """The annotated peaks from each sample are compared.
    3 Excel files (intensity.xlsx, mass.xlsx and mass_to_charge_and_intensity.xlsx) are generated.

    Args:
    peptide_file (string): Path to a FASTA-file with cyclotides of interest.
    """
    file_list = glob.glob('*.txt')
    file_list.sort(key=lambda x: int(x.split('.')[0][1::]))
    peak_cluster_list = list()

    fraction = 0
    for file in file_list:
        fraction += 1
        with open(file) as in_file:
            for line in in_file:
                line = line.split('\t')
                m_over_z = round(float(line[0]), 2)
                intensity = round(float(line[1]), 1)
                found = False
                for item in peak_cluster_list:
                    if fraction - item.max_fraction == 1:
                        if abs(item.peak_list[-1].my_m_over_z-m_over_z) <= 1:
                            item.peak_list.append(Peak(my_m_over_z=m_over_z, my_intensity=intensity, my_fraction=fraction))
                            found = True
                            item.max_fraction = fraction
                            break
                if not found:
                    peak_cluster_list.append(PeakCluster(max_fraction=fraction, peak_list=[Peak(my_m_over_z=m_over_z,
                                                                                                my_intensity=intensity,
                                                                                                my_fraction=fraction)]))

    peptide_dict = dict()
    peptides = SeqIO.parse(peptide_file, 'fasta')
    for peptide in peptides:
        name = peptide.name
        sequence = peptide.seq
        try:
            molecular_weight = 1-6*1.007825+Bio.SeqUtils.molecular_weight(sequence, seq_type='protein',
                                                                         circular=True, monoisotopic=True)
            peptide_dict[name] = molecular_weight
        except ValueError:
            print('Error in provided sequence file.')

    for cluster in peak_cluster_list:
        maximum_m_over_z = max(cluster.peak_list, key=lambda x: x.my_intensity).my_m_over_z
        for item in peptide_dict:
            if abs(peptide_dict[item]-maximum_m_over_z) <= 0.5:
                cluster.name_list.append(item)


    number_of_vertical_positions = len(peak_cluster_list)
    number_of_horizontal_positions = len(file_list)

    complete_lists_for_data_frame = list()
    intensity_lists_for_data_frame = list()
    mass_lists_for_data_frame = list()
    peak_cluster_list.sort()

    for cluster in peak_cluster_list:
        starting_value = cluster.peak_list[0].my_fraction
        end_value = cluster.peak_list[-1].my_fraction
        start = [' \n ']*(starting_value-1)
        end = [' \n ']*(number_of_horizontal_positions-end_value)
        complete_lists_for_data_frame.append(start+[str(i.my_m_over_z)+'\n'+str(i.my_intensity) for
                                                    i in cluster.peak_list]+end)
        start = [' ']*(starting_value-1)
        end = [' ']*(number_of_horizontal_positions-end_value)
        intensity_lists_for_data_frame.append(start+[i.my_intensity for i in cluster.peak_list]+end)
        mass_lists_for_data_frame.append(start+[i.my_m_over_z for i in cluster.peak_list]+end)


    index_list = [i.split('.')[0][1::] for i in file_list]
    peptide_name_list = list()
    for cluster in peak_cluster_list:
        peptide_name_list.append('\n'.join(cluster.name_list))


    df = pd.DataFrame(complete_lists_for_data_frame)
    df = df.transpose()
    df.index = index_list
    df.columns = peptide_name_list
    writer = pd.ExcelWriter('mass_to_charge_and_intensity.xlsx', engine='xlsxwriter')
    df.to_excel(writer, sheet_name='Sheet1', startrow=1, header=False)
    workbook = writer.book
    worksheet = writer.sheets['Sheet1']
    wrap_format = workbook.add_format({'text_wrap': True})
    for col_num, value in enumerate(df.columns.values):
        worksheet.write(0, col_num + 1, value, wrap_format)
    writer.save()


    df = pd.DataFrame(intensity_lists_for_data_frame)
    df = df.transpose()
    df.index = index_list
    df.columns = peptide_name_list
    writer = pd.ExcelWriter('intensity.xlsx', engine='xlsxwriter')
    df.to_excel(writer, sheet_name='Sheet1', startrow=1, header=False)
    workbook = writer.book
    worksheet = writer.sheets['Sheet1']
    wrap_format = workbook.add_format({'text_wrap': True})
    for col_num, value in enumerate(df.columns.values):
        worksheet.write(0, col_num + 1, value, wrap_format)
    writer.save()


    df = pd.DataFrame(mass_lists_for_data_frame)
    df = df.transpose()
    df.index = index_list
    df.columns = peptide_name_list
    writer = pd.ExcelWriter('mass.xlsx', engine='xlsxwriter')
    df.to_excel(writer, sheet_name='Sheet1', startrow=1, header=False)
    workbook = writer.book
    worksheet = writer.sheets['Sheet1']
    wrap_format = workbook.add_format({'text_wrap': True})
    for col_num, value in enumerate(df.columns.values):
        worksheet.write(0, col_num + 1, value, wrap_format)
    writer.save()



file_list = glob.glob('..\\..\\2021_03_31\\evaluated spectra\\*')
extract_peaks_from_msd_file(file_list=file_list)
compare_ms_spectra(peptide_file='peptides.fasta')
