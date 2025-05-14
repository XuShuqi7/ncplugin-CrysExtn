# NCrystal plugin CrysExtn

This plugin is dedicated to study the extinction effects observed in polycrystalline materials. Various models rooted in Darwin's mosaic crystal model picture are implemented, including Sabine's, Becker & Coppens', and Cooper & Rouse's extinction models. The random elastic deformation (RED) model developed by Kulda is also implemented. The syntax for using these models are presented as follows.

It is important here to clarify the meaning of crystallite and grain for the following parts. Under the mosaic crystal model picture of Darwin, the regions of perfect crystal lattices comprise a crystallite or perfect coherent domain. A grain of a polycrystalline material comprises an aggregate of crystallites. 

It should also be noted that extinction can be classified into two types: primary and secondary. Primary extinction arises from destructive interference between the incident and diffracted waves within a single crystallite, which reduces the diffracted beam intensity and consequently leads to an increase in transmission or a decrease in cross-section. Secondary extinction, on the other hand, involves multiple scattering among different crystallites with similar orientations. In this case, neutrons diffracted by one crystallite can be further diffracted by another, potentially redirecting them back toward the incident beam direction. Secondary extinction also accounts for the self-shielding effect among crystallites aligned along the neutron path, which reduces the overall diffracted intensity.

The extinction models implemented in this plugin have been investigated for the transmission measurements performed on four different beryllium grades using the HIPPO instrument at Los Alamos National Laboratory, and the experimental cross-sections were interpreted using Becker & Coppens' model. 

For more details, and if you use this plugin for your work, please refer and cite this publication: **Impact of extinction effects on neutron transmission and diffraction in solid beryllium metal, Journal of Applied Crystallography (under revision)**.

## Installation

Simply create a conda environment, install NCrystal via `conda install ncrystal`, and then `pip install "git+https://github.com/XuShuqi7/ncplugin-CrysExtn"` for installing the plugin

## Upgrade or reinstallation

`pip install --force-reinstall --upgrade --no-deps --no-build-isolation "git+https://github.com/XuShuqi7/ncplugin-CrysExtn"`

## Syntax

To use a plugin of NCrystal, one needs to include a supplymentrary section in the .ncmat file. This section is listed below for the aforementioned extinction models.

### Sabine's model

Sabine's model is further classified as uncorrelated and correlated block models. Here the term block is equivalent to crystallite. 

```
@CUSTOM_CRYSEXTN
  Sabine_uncorr  l  G  L  rect/tri
```
where `Sabine_uncorr` represents Sabine's uncorrelated block model which assumes no correlations between primary and secondary extinction effects, `l` is the crystallite size in unit of Aa, `G` represents the mosacity parameter in unit of 1/rad, L is the grain size (Aa), `rect/tri` represents the option for the distribution for the tilts between mosaic blocks, either rectangular or triangular distribution.

```
@CUSTOM_CRYSEXTN
  Sabine_corr  l  g  L
```
where `Sabine_corr` represents Sabine's correlated block model, `l` and `L` are the same physical quantities with same units as in the uncorrelated block model, `g` is also a mosacity parameter in unit of 1/rad.

### Becker & Coppens' (BC's) model

BC's extinction model is classified as both primary and secondary, pure primary, pure secondary, and pure secondary type-I and type-II.

```
@CUSTOM_CRYSEXTN
  BC_mix  l  g  L  Gauss/Lorentz/Fresnel
```
where `BC_mix` represents both primary and secondary extinction, `l` is the average path length through a crystallite (Aa), `g` is the mosacity parameter (1/rad), `L` is the average path length through a grain (Aa), and `Gauss/Lorentz/Fresnel` represents the option for calculating the quantities $A(\theta)$ and $B(\theta)$ involved in secondary extinction. For spherical domains, the average path length is related to the radius through the relations $l = (3/2)r$ and $L = (3/2)R$, for the crystallites and grains, respectively. In the paper, $l$ and $L$ are denoted as $t$ and $\overline{T}$, respectively.

```
@CUSTOM_CRYSEXTN
  BC_pure  l  g  L  Gauss/Lorentz/Fresnel
```
where `BC_pure` represents pure primary or pure type-I and type-II secondary extinction, `l`, `g`, `L`, and `Gauss/Lorentz/Fresnel` are the same as the primary and secondary extinction model. To trigger pure primary, one has to set `g=0` and `L=0`. For pure secondary type-I or type-II, one needs to set `l=0` and `g=0`, respectively.

```
@CUSTOM_CRYSEXTN
  BC_mod  l  g  L  Gauss/Lorentz/Fresnel
```
where `BC_mod` represents pure secondary extinction. The primary extinction factor $y_{\textrm{p}}$ has been fixed to 1, thus only secondary extinction is in operation. The other parameters are the same as the primary and secondary extinction model.

### Cooper & Rouse's (CR's) model

```
@CUSTOM_CRYSEXTN
  CR  l  g  L  
```
where `CR` represents Cooper & Rouse's extinction model, `l` is the crystallite size (Aa), g is the mosacity parameter (1/rad), and `L` is the grain size (Aa).


## References

### Sabine's model
  - International Tables for Crystallography (2006). Vol. C, Chapter 6.4, pp. 609–616.
### Becker & Coppens' model
  - Acta Cryst. (1974). A30, 129
  - Acta Cryst. (1995). A51, 662-667
### Cooper & Rouse's model
  - Acta Cryst. (1970). A26, 214
  - Acta Cryst. (1976). A32, 806

## Note

A jupyter notebook is attached in the folder "verification" for testing the implementation.
