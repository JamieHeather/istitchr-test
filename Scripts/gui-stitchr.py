import PySimpleGUI as sg
import functions as fxn
import stitchr as st
import collections as coll

# #create main window
# master = tkinter.Tk()
# master.title("tester")
# master.geometry("300x100")


# #make a label for the window
# label1 = tkinter.Label(master, text='Hellooooo')
# # Lay out label
# label1.pack()

# # Run forever!
# master.mainloop()

extra_gene_text = "(FASTA format, optional)"

col1 = [

    [sg.Button('Example data'), sg.Button('Reset form')],

    [sg.Text('Species', size=(30, 1), font=("Helvetica", 12))],
    [sg.Radio('Human', "RADIO1", key='rad_hs', default=True), sg.Radio('Mouse', "RADIO1", key='rad_mm')],

    [sg.Text('Additional genes', size=(30, 1), font=("Helvetica", 12))],
    [sg.MLine(default_text=extra_gene_text, size=(35, 3), key='additional_genes')],

    [sg.Button('Run Stitchr')],

    [sg.Button('Exit')]

]

col2 = [

    [sg.Text('TRA')],

    [sg.Text('TRAV')], [sg.InputText('', key='TRAV')],

    [sg.Text('TRAJ')], [sg.InputText('', key='TRAJ')],

    [sg.Text('TRA CDR3')], [sg.InputText('', key='TRA_cdr3')],

    [sg.Text('TRA out')], [sg.MLine(default_text='', size=(50, 20),
                                    key='TRA_out')]

]

col3 = [

    [sg.Text('TRB')],

    [sg.Text('TRBV')], [sg.InputText('', key='TRBV')],

    [sg.Text('TRBJ')], [sg.InputText('', key='TRBJ')],

    [sg.Text('TRB CDR3')], [sg.InputText('', key='TRB_cdr3')],

    [sg.Text('TRB out')], [sg.MLine(default_text='', size=(50, 20),
                                    key='TRB_out')]

]

layout = [[sg.Column(col1, element_justification='l'),
           sg.Column(col2, element_justification='l'),
           sg.Column(col3, element_justification='l')]]

window = sg.Window("stitchr", layout)
# event, values = window.read()

example_data = {
    'TRAV': 'TRAV1-2',
    'TRAJ': 'TRAJ33',
    'TRA_cdr3': 'CAVLDSNYQLIW',
    'TRA_out': '',
    'TRBV': 'TRBV7-3*01',
    'TRBJ': 'TRBJ1-1*01',
    'TRB_cdr3': 'CASSYLQAQYTEAFF',
    'TRB_out': '',
}

while True:
    event, values = window.read()

    if event == 'Example data':

        for field in example_data:
            window[field].update(example_data[field])

    elif event == 'Reset form':

        for field in example_data:
            window[field].update('')

    elif event == 'Run Stitchr':

        # Disable stitchr button while code is running
        window['Run Stitchr'].update(disabled=True)

        # Determine species
        if values['rad_hs']:
            species = 'HUMAN'
        elif values['rad_mm']:
            species = 'MOUSE'
        else:
            raise IOError('Species cannot be determined')

        print(species)

        # Loop through both chains, determine which are asked for, and read data in
        codons = fxn.get_optimal_codons('../Data/' + species + '/kazusa.txt', species)
        #        tcr_dat, functionality, out_lists, stitched, offsets, out_strs = [{}] * 6 # Todo rm
        outputs = coll.defaultdict()

        for chain in ['TRA', 'TRB']:

            if values[chain + 'V'] and values[chain + 'J'] and values[chain + '_cdr3']:
                tcr_dat, functionality = fxn.get_imgt_data(chain, st.gene_types, species)

                print(chain)
                # print(tcr_dat)

                tcr_bits = {'v': values[chain + 'V'], 'j': values[chain + 'J'], 'cdr3': values[chain + '_cdr3'],
                            'skip_c_checks': False, 'species': species, 'name': ''}
                tcr_bits = fxn.autofill_input(tcr_bits, chain)

                print(tcr_bits)

                # Run the stitching
                #                 out_lists[chain], stitched[chain], offsets[chain] = st.stitch(
                outputs[chain + '_out_list'], outputs[chain + '_stitched'], outputs[chain + '_offset'] = st.stitch(
                    tcr_bits, chain, tcr_dat, functionality, codons, 3)

                outputs[chain + '_out_str'] = '|'.join(outputs[chain + '_out_list'])
                outputs[chain + '_fasta'] = fxn.fastafy('nt|' + outputs[chain + '_out_str'],
                                                        outputs[chain + '_stitched'])

                window[chain + '_out'].update(outputs[chain + '_fasta'])

        #                 # Use the offset to 5' pad the stitched sequence with 'N's to make up for non-codon length 5' added sequences
        #                 print(fxn.fastafy('aa|' + out_str[chain], fxn.translate_nt('N' * offset[chain] + stitched[chain])))

        print('Stitching completed')

        # Re-enable stitchr button once completed
        window['Run Stitchr'].update(disabled=False)


    elif event in ('Exit', None):
        break

window.close()