# Sample Windows Indra Commands

You can literally copy these commands to your terminal to call Indra. 

**NB**: The first time, run either `RunThis_Windows.bat` or `RunThis_Unix.sh` on your machine (both located in the `installers` folder.

**NB**: I will use the command `python` everywhere, assuming that your default python installation is python3. If it isn't, or you are not sure, replace `python` with `python3` in all of the commands to be explicit and safe. 

## Train a model.

Using default options. At your terminal/command line use the command `python call_indra.py --help`.

### Windows
```
python3 call_indra.py --train 0 --stcode gen --n_samples 10 --fpath_in gen\gen_iwec.epw --fpath_out gen\gen_iwec_syn.epw --ftype epw --storepath gen

```

### Unix

```
python3 call_indra.py --train 0 --stcode gen --n_samples 10 --fpath_in gen/gen_iwec.epw --fpath_out gen\gen_iwec_syn.epw --ftype epw --storepath gen

```
