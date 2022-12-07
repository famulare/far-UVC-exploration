# Repo for stuff that's helping me learn about far UV-C sterilization.

## My twitter threads
- 6 December 2022: [Thread tree starting from Siyuan Wang's "spherical chicken" 3D room modeling](https://mobile.twitter.com/SiyuanWang_PhD/status/1598966798334201856). Follow all the branching to see the full discussion. 
    - I want to highlight [a discussion of jet dynamics from breathing and talking](https://mobile.twitter.com/famulare_mike/status/1600262799103758336) informed by a fantastic review, [Effects of ventilation on the indoor spread of COVID-19 by Bhagat, Davies Wykes, Dalziel, and Linden](https://www.cambridge.org/core/journals/journal-of-fluid-mechanics/article/effects-of-ventilation-on-the-indoor-spread-of-covid19/CF272DAD7C27DC44F6A9393B0519CAE3) and its [excellent videos](https://www.cambridge.org/core/journals/journal-of-fluid-mechanics/article/effects-of-ventilation-on-the-indoor-spread-of-covid19/CF272DAD7C27DC44F6A9393B0519CAE3#supplementary-materials). 
    - The conclusion we're heading toward is that trying to create a personal bubble with low power, easily portable far UV-C in a dynamic environment is unlikely to be effective. (Would like to be wrong though...)
    - However, that does not negate the value of airborne mitigations to reduce far-field transmission, as [reiterated here](https://mobile.twitter.com/famulare_mike/status/1600282668490690560).
- 7 November 2022: [I threw together a crude, amateur model of the beam cone for an Ergo X One personal far-UV-C device. With that, a few thoughts...](https://twitter.com/famulare_mike/status/1589775727754637313)
- 2 October 2022: [Very amateur notes on the UV Can Lily personal far UV-C device](https://twitter.com/famulare_mike/status/1576736411616772096)
- 23 September 2022: Early pessimism based on sloppy 1/r^2 reasoning. Has good questions but disregard take on the Ergo X One. ["I love this idea but I’m not seeing the math work out...."](https://twitter.com/famulare_mike/status/1573514122075181056)

## Devices considered
- [UV Can Lily](https://uv-can.com/products/lily-handheld-personal-far-uv-disinfection-light)
- [UV Medico UV222](https://uvmedico.com/lamps/uv222-specifications/)
- [Ergo X One](https://www.ergo-healthtech.com/xone)

## What's up with `DEPRECATED`?
All files with `CRUDE_SKETCH` in the name are based on digital spectrometer or manufacturer-reported irradiance data. The `DEPRECATED_ergo_x_one_ushio_222_card_experiment` directory for the Ergo X One is based on estimating irradiance using [Ushio Dose222 Far UV-C Indicator cards](https://www.ushio.com/products/infection-prevention/dose222-far-uv-c-indicator-cards/##). It was impossible to match colors by eye with precision, and compared to digital spectrometer readings, these irradiance estimates are 2-3x low, and so should be ignored. I'm including this in the repo in case I want to repeat or improve this experiment, which is much cheaper than buying a spectrometer if you only want to do it once.
