import primer3
import random
import os

values_to_write = [
    'SEQUENCE_PRIMER',
    'PRIMER_LEFT_0_TM',
    'PRIMER_LEFT_0_GC_PERCENT',
    'PRIMER_LEFT_0_SELF_ANY_TH',
    'PRIMER_LEFT_0_SELF_END_TH',
    'PRIMER_LEFT_0_HAIRPIN_TH',
    'PRIMER_RUNS',
    'PRIMER_DIMER_OK',
    'PRIMER_PARAMS_OK'
]
primer3_path = 'primer3_core'


class PrimerGenerator:
    def __init__(self, check_runs: bool, check_tm: bool, check_gc: bool, tm_max: float, tm_min: float, gc_max: int, gc_min: int) -> None:
        self.check_runs = check_runs
        self.check_tm = check_tm
        self.check_gc = check_gc
        if self.check_gc:
            self.gc_min = gc_min
            self.gc_max = gc_max
        if self.check_tm:
            self.tm_min = tm_min
            self.tm_max = tm_max
        pass


    def primer_check_runs(self, primer: str) -> bool:
        return not ('AAAAA' in primer or 'CCCCC' in primer or 'TTTTT' in primer or 'GGGGG' in primer)


    def primer_check_temp(self, primer: str) -> bool:
        tm = primer3.calcTm(primer, mv_conc=50, dv_conc=4, dntp_conc=0.5, dna_conc=50)
        return self.tm_min < tm < self.tm_max


    def primer_check_gc(self, primer: str) -> bool:
        gc = ((primer.count('C') + primer.count('G')) / float(len(primer)) )* 100
        return self.gc_min < gc < self.gc_max


    def generate_primers(self, count: int, lenght: int):
        primers = []
        generated = 0
        while generated < count:
            primer = ''.join(random.choice('ACTG') for _ in range(lenght))
            primers.append(primer)
            generated += 1
        return primers

if __name__ == '__main__':
    primer_count = int(input('Ile primerów wygenerować? '))
    primer_len = int(input('Długość primera? '))
    # primer_check_runs_input = input('Sprawdzać runs? [T/n] ')
    # if primer_check_runs_input.lower().strip() == 't':
    #     primer_check_runs = True
    # else:
    #     primer_check_runs = False
    # primer_check_temp_input = input('Sprawdzać temperaturę topnienia? [T/n] ')
    # if primer_check_temp_input.lower().strip() == 't':
    #     primer_check_temp = True
    #     primer_tm_min = float(input('\tMinimalna temperatura topnienia? '))
    #     primer_tm_max = float(input('\tMaksymalna temperatura topnienia? '))
    # else:
    #     primer_check_temp = False
    #     primer_tm_min = None
    #     primer_tm_max = None
    # primer_check_gc_input = input('Sprawdzać zawartość GC? [T/n] ')
    # if primer_check_gc_input.lower().strip() == 't':
    #     primer_check_gc = True
    #     primer_gc_min = float(input('\tMinimalna zawartość GC? '))
    #     primer_gc_max = float(input('\tMaksymalna zawartość GC? '))
    # else:
    #     primer_check_gc = False
    #     primer_gc_min = None
    #     primer_gc_max = None
    output_file_name = input('Nazwa pliku wyjściowego? ')

    print('Generowanie primerów...')
    # primer_generator = PrimerGenerator(primer_check_runs, primer_check_temp, primer_check_gc, primer_tm_max, primer_tm_min, primer_gc_max, primer_gc_min)
    primer_generator = PrimerGenerator(False, False, False, None, None, None, None)
    primers = primer_generator.generate_primers(primer_count, primer_len)

    print('\nPrzygotowywanie danych wejściowych dla primer3_core...')
    with open('p3_input', 'w') as p3_input_file:
        psm = 50.0
        psd = 4.0
        pdc = 0.5
        pnc = 50.0
        end = '\n'
        p3_input_file.writelines([
            'SEQUENCE_ID=temp_file' + end,
            'PRIMER_TASK=check_primers' + end,
            'PRIMER_SALT_MONOVALENT=' + str(psm) + end,
            'PRIMER_SALT_DIVALENT=' + str(psd) + end,
            'PRIMER_DNTP_CONC=' + str(pdc) + end,
            'PRIMER_DNA_CONC=' + str(pnc) + end
        ])
        for primer in primers:
                p3_input_file.write('SEQUENCE_PRIMER=' + primer + end)
                p3_input_file.write('=\n')
        p3_input_file.close()

    print('Analiza primerów przez primer3_core...')
    os.system(primer3_path + ' --output=p3_output < p3_input')

    print('Przetwarzanie danych wyjściowych primer3_core na csv...')
    primers_csv_data = ''
    generated_primers = 0
    ok_primers = 0
    with open('p3_output', 'r') as p3_output_file:
        primers_data = p3_output_file.read().split('=\n')
        for primer_data in primers_data:
            primer_dict = {}
            for line in primer_data.split('\n'):
                if len(line) and len(line.split('=')) == 2:
                    [key, value] = line.split('=')
                    primer_dict[str(key).strip()] = value.strip()
            if primer_dict.keys():
                try:
                    primer_dimer_ok = '1' if sum(float(primer_dict[k]) for k in values_to_write[-6:-3]) < 0.5 else '0'
                    primer_params_ok = '1' if 52.0 < float(primer_dict['PRIMER_LEFT_0_TM']) < 58.0 and 45.0 < float(primer_dict['PRIMER_LEFT_0_GC_PERCENT']) < 60.0 else '0'
                    primer_runs = sum(primer_dict['SEQUENCE_PRIMER'].count(run) for run in ['AAAA', 'CCCC', 'TTTT', 'GGGG'])
                    primer_dict['PRIMER_DIMER_OK'] = primer_dimer_ok
                    primer_dict['PRIMER_PARAMS_OK'] = primer_params_ok
                    primer_dict['PRIMER_RUNS'] = str(primer_runs)
                    generated_primers = generated_primers + 1
                    primers_csv_data = primers_csv_data + ','.join(primer_dict[k] for k in values_to_write) + end
                    if primer_dimer_ok == '1' and primer_params_ok == '1' and primer_runs == 0:
                        ok_primers += 1
                except KeyError:
                    pass
        p3_output_file.close()
    
    print('Zapisywanie danych do pliku wyjściowego...')
    with open(output_file_name, 'w') as output_file:
        output_file.write(','.join(values_to_write) + end)
        output_file.write(primers_csv_data)
        output_file.close()
    # os.remove('p3_input')
    # os.remove('p3_output')
    print('Zakończono. Wygenerowano %d primerów - %d oznaczonych jako OK.' % (generated_primers, ok_primers))




