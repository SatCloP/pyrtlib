
from pyrtlib import __version__
from pyrtlib.absorption_model import H2OAbsModel
from pyrtlib.atmospheric_profiles import AtmosphericProfiles as atmp
from pyrtlib.tb_spectrum import TbCloudRTE
from pyrtlib.utils import ppmv2gkg, mr2rh, get_frequencies_sat


def main():
    z, p, _, t, md = atmp.gl_atm(atmp.TROPICAL)

    gkg = ppmv2gkg(md[:, atmp.H2O], atmp.H2O)
    rh = mr2rh(p, t, gkg)[0] / 100

    frq = get_frequencies_sat("MWI")

    rte = TbCloudRTE(z, p, t, rh, frq)
    rte.init_absmdl('R20')
    H2OAbsModel.model = 'R21SD'
    H2OAbsModel.set_ll()
    df = rte.execute()

    print('Hello Spectrum!')
    print(df)
    print(f'PyRTLib {__version__} successfully installed')


if __name__ == '__main__':
    main()
