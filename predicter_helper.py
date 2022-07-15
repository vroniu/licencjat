from keras.models import load_model, Sequential
from converter import convert_x, convert_y_predicted
import pandas as pd


def load_cnn_model():
    model = load_model('cnn_model')
    model.load_weights('cnn_model_weights')
    return model


def convert_seq_to_model_input(sequence):
    translate_dict = {
        ord("A"):ord("1"),
        ord("C"):ord("2"),
        ord("T"):ord("3"),
        ord("G"):ord("4")
    }
    sequence = sequence.translate(translate_dict)
    dataframe = pd.DataFrame(list(sequence)).T
    return convert_x(dataframe)

def prepare_primer3_input(template, target_start, amplified_seq_len):
    input_string = """PRIMER_SALT_MONOVALENT=50.0
PRIMER_SALT_DIVALENT=4.0
PRIMER_DNTP_CONC=0.5
PRIMER_DNA_CONC=50.0
PRIMER_TASK=generic
SEQUENCE_TEMPLATE=ENTER_SEQUENCE_HERE
SEQUENCE_TARGET=ENTER_TARGET_START,ENTER_AMPLIFIED_SEQ_LEN
="""
    return input_string.replace("ENTER_SEQUENCE_HERE", template).replace("ENTER_TARGET_START", str(target_start)).replace("ENTER_AMPLIFIED_SEQ_LEN", str(amplified_seq_len))


def get_primer_parameters(lines, primer_index):
    left_primer = next(line.replace("PRIMER_LEFT_%d_SEQUENCE=" % primer_index, "") for line in lines if line.startswith("PRIMER_LEFT_%d_SEQUENCE=" % primer_index))
    right_primer = next(line.replace("PRIMER_RIGHT_%d_SEQUENCE=" % primer_index, "") for line in lines if line.startswith("PRIMER_RIGHT_%d_SEQUENCE=" % primer_index))
    left_primer_start = next(line.replace("PRIMER_LEFT_%d=" % primer_index, "") for line in lines if line.startswith("PRIMER_LEFT_%d=" % primer_index))
    right_primer_start = next(line.replace("PRIMER_RIGHT_%d=" % primer_index, "") for line in lines if line.startswith("PRIMER_RIGHT_%d=" % primer_index))
    return (left_primer_start, left_primer, right_primer_start, right_primer)

def parse_primer3_output(raw_output):
    result = []
    lines = raw_output.split("\n")
    primers_returned = int(next(line.replace("PRIMER_PAIR_NUM_RETURNED=", "") for line in lines if line.startswith("PRIMER_PAIR_NUM_RETURNED")))
    for i in range(primers_returned):
        result.append(get_primer_parameters(lines, i))
    return result

