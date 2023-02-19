This is the repository for CLM-Palm development initiated by Yuanchao Fan based on the following references:

Fan, Y., 2016. Modeling oil palm monoculture and its associated impacts on land-atmosphere carbon, water and energy fluxes in Indonesia (Doctoral thesis). University of Göttingen, Germany.

Fan, Y., Roupsard, O., Bernoux, M., Le Maire, G., Panferov, O., Kotowska, M.M., Knohl, A., 2015. A sub-canopy structure for simulating oil palm in the Community Land Model (CLM-Palm): phenology, allocation and yield. Geoscientific Model Development 8, 3785–3800.

Fan, Y., Meijide, A., Lawrence, D.M., Roupsard, O., Carlson, K.M., Chen, H.-Y., Röll, A., Niu, F., Knohl, A., 2019. Reconciling Canopy Interception Parameterization and Rainfall Forcing Frequency in the Community Land Model for Simulating Evapotranspiration of Rainforests and Oil Palm Plantations in Indonesia. Journal of Advances in Modeling Earth Systems 11, 732–751.

The repository also includes code of a perennial semi-deciduous woody crop PFT for rubber, described in the following paper:

Ali, A.A., Fan, Y., Corre, M.D., Kotowska, M.M., Preuss-Hassler, E., Cahyo, A.N., Moyano, F.E., Stiegler, C., Röll, A., Meijide, A., Olchev, A., Ringeler, A., Leuschner, C., Ariani, R., June, T., Tarigan, S., Kreft, H., Hölscher, D., Xu, C., Koven, C.D., Dagon, K., Fisher, R.A., Veldkamp, E., Knohl, A., 2022. Implementing a New Rubber Plant Functional Type in the Community Land Model (CLM5) Improves Accuracy of Carbon and Water Flux Estimation. Land 11, 183.

The CLM-Palm repository has undergone upgrades from svn to git and from the base CLM4.5 to CLM5.0. The active branch "clm-palm-dev01-ctsm1.0.dev040" is the latest version based on the trunk of ctsm1.0-dev (similar to the trunk clm5.0 release).

Main features of CLM-Palm include a perennial crop phenology and carbon/nitrogen allocation for oil palm plantations.

There are two options of perennial phenology and canopy structure for oil palm:

1. A generic perennial evergreen crop phenology that represents consecutive fruit filling and fruit harvests at monthly or yearly intervals depending on parameterization.

This perennial phenology not only works for oil palm but also for other perennial crops such as rubber, coca, coconut, coffee etc.

The generic perennial phenology code is embedded in the CropPhenology subroutine in CNPhenologyMod, and is activated by setting perennial=1 and phytomer=0.

2. A multilayer canopy phytomer-based oil palm structure and phenology and allocation subroutines as described in Fan et al. (2015). This phytomer-based version allows to simulate oil palm's growth, yield, and management (leaf pruning and fruit harvest) more realistically according to observations.

The phenology code is in an independent subroutine called PalmPhenology in CNPhenologyMod and is activated by setting perennial=1 and phytomer=1.

Phytomer specific leaf and fruit allocations are in clauses where phytomer>0 or mxnp>0 in NutrientCompetitionCLM45defaultMod.F90 and NutrientCompetitionFlexibleCNMod.F90

A mutlilayer radiative transfer code is also developed for oil palm and applicable to other tree PFTs, which is described in Fan 2016 and not published in journals yet.

For collaborations and questions on this repository, please write to Yuanchao (yuanchao.fan@sz.tsinghua.edu.cn, yfansunny@gmail.com), or to Ashehad (ashehad.ali@uni-goettingen.de)




Documentation of CTSM:

The Community Terrestrial Systems Model.
This includes the Community Land Model (CLM5.0 and CLM4.5) of the Community Earth System Model.
For documentation, quick start, diagnostics, model output and
references, see
http://www.cesm.ucar.edu/models/cesm2.0/land/
and
https://escomp.github.io/ctsm-docs/
For help with how to work with CTSM in git, see
https://github.com/ESCOMP/ctsm/wiki/Getting-started-with-CTSM-in-git
and
https://github.com/ESCOMP/ctsm/wiki/Recommended-git-setup
To get updates on CTSM tags and important notes on CTSM developments
join our low traffic email list:
https://groups.google.com/a/ucar.edu/forum/#!forum/ctsm-dev
(Send email to ctsm-software@ucar.edu if you have problems with any of this)
