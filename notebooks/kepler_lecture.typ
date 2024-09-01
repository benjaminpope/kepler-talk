// Some definitions presupposed by pandoc's typst output.
#let blockquote(body) = [
  #set text( size: 0.92em )
  #block(inset: (left: 1.5em, top: 0.2em, bottom: 0.2em))[#body]
]

#let horizontalrule = [
  #line(start: (25%,0%), end: (75%,0%))
]

#let endnote(num, contents) = [
  #stack(dir: ltr, spacing: 3pt, super[#num], contents)
]

#show terms: it => {
  it.children
    .map(child => [
      #strong[#child.term]
      #block(inset: (left: 1.5em, top: -0.4em))[#child.description]
      ])
    .join()
}

// Some quarto-specific definitions.

#show raw.where(block: true): block.with(
    fill: luma(230), 
    width: 100%, 
    inset: 8pt, 
    radius: 2pt
  )

#let block_with_new_content(old_block, new_content) = {
  let d = (:)
  let fields = old_block.fields()
  fields.remove("body")
  if fields.at("below", default: none) != none {
    // TODO: this is a hack because below is a "synthesized element"
    // according to the experts in the typst discord...
    fields.below = fields.below.amount
  }
  return block.with(..fields)(new_content)
}

#let empty(v) = {
  if type(v) == "string" {
    // two dollar signs here because we're technically inside
    // a Pandoc template :grimace:
    v.matches(regex("^\\s*$")).at(0, default: none) != none
  } else if type(v) == "content" {
    if v.at("text", default: none) != none {
      return empty(v.text)
    }
    for child in v.at("children", default: ()) {
      if not empty(child) {
        return false
      }
    }
    return true
  }

}

// Subfloats
// This is a technique that we adapted from https://github.com/tingerrr/subpar/
#let quartosubfloatcounter = counter("quartosubfloatcounter")

#let quarto_super(
  kind: str,
  caption: none,
  label: none,
  supplement: str,
  position: none,
  subrefnumbering: "1a",
  subcapnumbering: "(a)",
  body,
) = {
  context {
    let figcounter = counter(figure.where(kind: kind))
    let n-super = figcounter.get().first() + 1
    set figure.caption(position: position)
    [#figure(
      kind: kind,
      supplement: supplement,
      caption: caption,
      {
        show figure.where(kind: kind): set figure(numbering: _ => numbering(subrefnumbering, n-super, quartosubfloatcounter.get().first() + 1))
        show figure.where(kind: kind): set figure.caption(position: position)

        show figure: it => {
          let num = numbering(subcapnumbering, n-super, quartosubfloatcounter.get().first() + 1)
          show figure.caption: it => {
            num.slice(2) // I don't understand why the numbering contains output that it really shouldn't, but this fixes it shrug?
            [ ]
            it.body
          }

          quartosubfloatcounter.step()
          it
          counter(figure.where(kind: it.kind)).update(n => n - 1)
        }

        quartosubfloatcounter.update(0)
        body
      }
    )#label]
  }
}

// callout rendering
// this is a figure show rule because callouts are crossreferenceable
#show figure: it => {
  if type(it.kind) != "string" {
    return it
  }
  let kind_match = it.kind.matches(regex("^quarto-callout-(.*)")).at(0, default: none)
  if kind_match == none {
    return it
  }
  let kind = kind_match.captures.at(0, default: "other")
  kind = upper(kind.first()) + kind.slice(1)
  // now we pull apart the callout and reassemble it with the crossref name and counter

  // when we cleanup pandoc's emitted code to avoid spaces this will have to change
  let old_callout = it.body.children.at(1).body.children.at(1)
  let old_title_block = old_callout.body.children.at(0)
  let old_title = old_title_block.body.body.children.at(2)

  // TODO use custom separator if available
  let new_title = if empty(old_title) {
    [#kind #it.counter.display()]
  } else {
    [#kind #it.counter.display(): #old_title]
  }

  let new_title_block = block_with_new_content(
    old_title_block, 
    block_with_new_content(
      old_title_block.body, 
      old_title_block.body.body.children.at(0) +
      old_title_block.body.body.children.at(1) +
      new_title))

  block_with_new_content(old_callout,
    new_title_block +
    old_callout.body.children.at(1))
}

// 2023-10-09: #fa-icon("fa-info") is not working, so we'll eval "#fa-info()" instead
#let callout(body: [], title: "Callout", background_color: rgb("#dddddd"), icon: none, icon_color: black) = {
  block(
    breakable: false, 
    fill: background_color, 
    stroke: (paint: icon_color, thickness: 0.5pt, cap: "round"), 
    width: 100%, 
    radius: 2pt,
    block(
      inset: 1pt,
      width: 100%, 
      below: 0pt, 
      block(
        fill: background_color, 
        width: 100%, 
        inset: 8pt)[#text(icon_color, weight: 900)[#icon] #title]) +
      if(body != []){
        block(
          inset: 1pt, 
          width: 100%, 
          block(fill: white, width: 100%, inset: 8pt, body))
      }
    )
}



#let article(
  title: none,
  authors: none,
  date: none,
  abstract: none,
  cols: 1,
  margin: (x: 1.25in, y: 1.25in),
  paper: "us-letter",
  lang: "en",
  region: "US",
  font: (),
  fontsize: 11pt,
  sectionnumbering: none,
  toc: false,
  toc_title: none,
  toc_depth: none,
  toc_indent: 1.5em,
  doc,
) = {
  set page(
    paper: paper,
    margin: margin,
    numbering: "1",
  )
  set par(justify: true)
  set text(lang: lang,
           region: region,
           font: font,
           size: fontsize)
  set heading(numbering: sectionnumbering)

  if title != none {
    align(center)[#block(inset: 2em)[
      #text(weight: "bold", size: 1.5em)[#title]
    ]]
  }

  if authors != none {
    let count = authors.len()
    let ncols = calc.min(count, 3)
    grid(
      columns: (1fr,) * ncols,
      row-gutter: 1.5em,
      ..authors.map(author =>
          align(center)[
            #author.name \
            #author.affiliation \
            #author.email
          ]
      )
    )
  }

  if date != none {
    align(center)[#block(inset: 1em)[
      #date
    ]]
  }

  if abstract != none {
    block(inset: 2em)[
    #text(weight: "semibold")[Abstract] #h(1em) #abstract
    ]
  }

  if toc {
    let title = if toc_title == none {
      auto
    } else {
      toc_title
    }
    block(above: 0em, below: 2em)[
    #outline(
      title: toc_title,
      depth: toc_depth,
      indent: toc_indent
    );
    ]
  }

  if cols == 1 {
    doc
  } else {
    columns(cols, doc)
  }
}

#set table(
  inset: 6pt,
  stroke: none
)
#show: doc => article(
  title: [Interactive Data Analysis with Kepler/K2 and TESS],
  toc_title: [Table of contents],
  toc_depth: 3,
  cols: 1,
  doc,
)


#strong[Benjamin Pope]

#emph[University of Queensland]

#block[
```python
# This is going to be an interactive lecture using Jupyter! 
# So let's import things!

import numpy as np 
import matplotlib.pyplot as plt

import lightkurve
import k2sc

from astropy.timeseries import LombScargle, BoxLeastSquares
import astropy.units as u

import warnings; warnings.simplefilter('ignore')

%matplotlib inline
```

]
= Kepler data
<kepler-data>
The Kepler Space Telescope revolutionised exoplanetary science by providing high-precision photometry of hundreds of thousands of stars. The data is publicly available and has been used to discover thousands of exoplanets.

While the mission science team originally aimed to keep their data proprietary and analyse it in-house, after a brief period they pivoted to an open-source model with data hosted on (MAST)\[https:\/\/archive.stsci.edu/kepler/\], and many teams developed open-source software to analyse the data. This is what I did my DPhil on!

We’re going to interactively explore data from Kepler to see how we would discover planets, stellar variability, and correct for instrumental systematics.

```python
# download Kepler-10 data
target = 'Kepler-10'
search = lightkurve.search_targetpixelfile(target,exptime=60) # get 1 minute cadence data
search 
```

```python
tpf = search[1].download()
tpf.plot()
```

#box(image("kepler_lecture_files/figure-typst/cell-4-output-1.png"))

Cool! We can make this into a light curve the simplest way possible: define a pixel mask, and sum the flux in the pixels.

Let’s use the default mask from the pipeline.

```python
tpf.plot(aperture_mask = 'pipeline')
```

#box(image("kepler_lecture_files/figure-typst/cell-5-output-1.png"))

```python
# let's sum this to make a lightcurve
lc = tpf.to_lightcurve(aperture_mask = 'pipeline').normalize().remove_nans()
lc = lc[lc.quality==0]
lc.scatter()
```

#box(image("kepler_lecture_files/figure-typst/cell-6-output-1.png"))

```python
# what about that long-term trend?
flat = lc.flatten(window_length=2*24*60+1).remove_outliers(sigma_upper=3,sigma_lower=100)
flat.plot()
```

#box(image("kepler_lecture_files/figure-typst/cell-7-output-1.png"))

```python
our_lc = flat.bin(0.25*u.hour).remove_nans()
our_lc.plot()
```

#box(image("kepler_lecture_files/figure-typst/cell-8-output-1.png"))

#block[
```python
# box-least-squares fits a periodic signal of dips to the light curve

bls = BoxLeastSquares(our_lc.time, our_lc.flux, dy=our_lc.flux_err)

periods = np.linspace(0.3, 1.5, 1000)*u.day
durations = np.linspace(0.05, 0.2, 10)*u.day
periodogram = bls.power(periods, durations)

best_period = periods[np.argmax(periodogram.depth_snr.value)]
best_duration = periodogram.duration[np.argmax(periodogram.depth_snr.value)]
best_depth = periodogram.depth[np.argmax(periodogram.depth_snr.value)]
t0 = periodogram.transit_time[np.argmax(periodogram.depth_snr.value)]
```

]
```python
periodogram.keys()
```

```
dict_keys(['objective', 'period', 'power', 'depth', 'depth_err', 'duration', 'transit_time', 'depth_snr', 'log_likelihood'])
```

```python
plt.plot(periodogram.period.value,periodogram.depth_snr.value)
plt.xlabel('Period [days]')
plt.ylabel('Signal-to-noise')
plt.xlim(periods.min().value, periods.max().value)
plt.axvline(best_period.value, color='red', linestyle='dashed')
plt.text(best_period.value+0.05, np.max(periodogram.depth_snr.value), f'Best period: {best_period.value:.4g} days')
```

```
Text(0.8869369369369369, 40.83467531875844, 'Best period: 0.8369 days')
```

#block[
#box(image("kepler_lecture_files/figure-typst/cell-11-output-2.png"))

]
```python
# what does that model look like?
model = bls.model(our_lc.time, best_period, best_duration, t0)
plt.plot(our_lc.time.value, our_lc.flux.value)
plt.plot(our_lc.time.value, model, color='red', lw=1)
```

#box(image("kepler_lecture_files/figure-typst/cell-12-output-1.png"))

```python
# fold it
folded_lc = our_lc.fold(period=best_period, t0=t0)

folded_times = our_lc.time.value % best_period.value
folded_args = np.argsort(folded_times)

plt.plot(folded_times[folded_args],our_lc.flux.value[folded_args]-1,'.k')
plt.plot(folded_times[folded_args],model[folded_args]-1,color='C3',lw=4)
plt.xlim(0, best_period.value)
plt.text(0.35, -best_depth*2, f'Best Depth: {best_depth*1e6:.0f} ppm',color='C3',fontdict={'size':14,'weight':'bold'})
```

```
Text(0.35, -0.00031383340219341395, 'Best Depth: 157 ppm')
```

#block[
#box(image("kepler_lecture_files/figure-typst/cell-13-output-2.png"))

]
#block[
```python
# what is Rp/Rs?
print(f'Rp/R* =  {(best_depth.value)**0.5 * 100:.2f} %')
# Earth for comparison is 0.00916794
print(f'Earth/Sun Radius =  {0.00916794 * 100:.2f} %')
```

#block[
```
Rp/R* =  1.25 %
Earth/Sun Radius =  0.92 %
```

]
]
Now let’s look at how it was in K2: where the systematics are terrible! Let’s look at my favourite star, EPIC 212521166.

```python
k2target = 'EPIC 212521166'
search = lightkurve.search_lightcurvefile(k2target)
search
```

```python
lc = search.download().PDCSAP_FLUX.normalize().remove_nans()
lc = lc[lc.quality==0]

lc.scatter()
```

#box(image("kepler_lecture_files/figure-typst/cell-16-output-1.png"))

Look at all those transits! But also look how it has a long-term variation, which is probably stellar, and a short term variation from the 6-hour pointing jitter of K2. We’re going to correct for this by #emph[detrending] with a Gaussian process.

```python
x, y = lc.pos_corr1, lc.pos_corr2

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4))

ax1.plot(lc.time.value, x, 'k.', markersize=3)
ax1.plot(lc.time.value, y, 'r.', markersize=3)
ax1.set_xlabel('Time (BJD)')
ax1.set_ylabel('Position')
ax1.legend(['X', 'Y'])

ax2.plot(x, y, '.k', markersize=3)
ax2.set_xlabel('X')
ax2.set_ylabel('Y')
```

```
Text(0, 0.5, 'Y')
```

#block[
#box(image("kepler_lecture_files/figure-typst/cell-17-output-2.png"))

]
#block[
```python
from k2sc.standalone import k2sc_lc


lc.__class__ = k2sc_lc
lc.k2sc()
```

#block[
```
/Users/benpope/opt/anaconda3/envs/lk/lib/python3.10/site-packages/k2sc/standalone.py:293: LightkurveDeprecationWarning: The hdu function is deprecated and may be removed in a future version.
        Use fits.open(lc.filename) instead.
  primary_header = self.hdu[0].header
/Users/benpope/opt/anaconda3/envs/lk/lib/python3.10/site-packages/k2sc/standalone.py:293: ResourceWarning: unclosed file <_io.BufferedReader name='/Users/benpope/.lightkurve-cache/mastDownload/K2/ktwo212521166-c06_lc/ktwo212521166-c06_llc.fits'>
  primary_header = self.hdu[0].header
/Users/benpope/opt/anaconda3/envs/lk/lib/python3.10/site-packages/k2sc/standalone.py:294: LightkurveDeprecationWarning: The hdu function is deprecated and may be removed in a future version.
        Use fits.open(lc.filename) instead.
  data_header = self.hdu[1].header
/Users/benpope/opt/anaconda3/envs/lk/lib/python3.10/site-packages/k2sc/standalone.py:294: ResourceWarning: unclosed file <_io.BufferedReader name='/Users/benpope/.lightkurve-cache/mastDownload/K2/ktwo212521166-c06_lc/ktwo212521166-c06_llc.fits'>
  data_header = self.hdu[1].header
/Users/benpope/opt/anaconda3/envs/lk/lib/python3.10/site-packages/k2sc/dtdata.py:63: DeprecationWarning: `np.int` is a deprecated alias for the builtin `int`. To silence this warning, use `int` by itself. Doing this will not modify any behavior and is safe. When replacing `np.int`, you may wish to use e.g. `np.int64` or `np.int32` to specify the precision. If you wish to review your current use, check the release note link for additional information.
Deprecated in NumPy 1.20; for more details and guidance: https://numpy.org/devdocs/release/1.20.0-notes.html#deprecations
  bstarts = np.full(nblocks, -bspan, np.int)   ## Starting indices for blocks
/Users/benpope/opt/anaconda3/envs/lk/lib/python3.10/site-packages/k2sc/dtdata.py:17: DeprecationWarning: `np.bool` is a deprecated alias for the builtin `bool`. To silence this warning, use `bool` by itself. Doing this will not modify any behavior and is safe. If you specifically wanted the numpy scalar type, use `np.bool_` here.
Deprecated in NumPy 1.20; for more details and guidance: https://numpy.org/devdocs/release/1.20.0-notes.html#deprecations
  self._mask   = array(mask) if mask is not None else ones(self._flux.size, np.bool)
/Users/benpope/opt/anaconda3/envs/lk/lib/python3.10/site-packages/k2sc/gp.py:98: DeprecationWarning: `np.bool` is a deprecated alias for the builtin `bool`. To silence this warning, use `bool` by itself. Doing this will not modify any behavior and is safe. If you specifically wanted the numpy scalar type, use `np.bool_` here.
Deprecated in NumPy 1.20; for more details and guidance: https://numpy.org/devdocs/release/1.20.0-notes.html#deprecations
  mask = zeros((t1.size,t2.size), np.bool)
```

]
#block[
```
Using default splits [2390, 2428] for campaign 6
Starting initial outlier detection
  Flagged 34 ( 1.0%) outliers.
Starting Lomb-Scargle period search
  Using SqrExp position kernel
  Found periodicity p =   14.79 (fap 1.7550e-190 < 1e-50), will use a quasiperiodic kernel
Starting global hyperparameter optimisation using DE
```

]
#block[
```
  0%|          | 0/150 [00:00<?, ?it/s]/Users/benpope/opt/anaconda3/envs/lk/lib/python3.10/site-packages/k2sc/de.py:85: DeprecationWarning: `np.int` is a deprecated alias for the builtin `int`. To silence this warning, use `int` by itself. Doing this will not modify any behavior and is safe. When replacing `np.int`, you may wish to use e.g. `np.int64` or `np.int32` to specify the precision. If you wish to review your current use, check the release note link for additional information.
Deprecated in NumPy 1.20; for more details and guidance: https://numpy.org/devdocs/release/1.20.0-notes.html#deprecations
  t = np.zeros(3, np.int)
 68%|██████▊   | 102/150 [04:57<02:19,  2.91s/it, -ln(L)=-6018.5]
/Users/benpope/opt/anaconda3/envs/lk/lib/python3.10/site-packages/k2sc/gp.py:98: DeprecationWarning: `np.bool` is a deprecated alias for the builtin `bool`. To silence this warning, use `bool` by itself. Doing this will not modify any behavior and is safe. If you specifically wanted the numpy scalar type, use `np.bool_` here.
Deprecated in NumPy 1.20; for more details and guidance: https://numpy.org/devdocs/release/1.20.0-notes.html#deprecations
  mask = zeros((t1.size,t2.size), np.bool)
```

]
#block[
```
  DE finished in 297 seconds
  DE minimum found at: [-5.997e+00  9.718e-01  1.254e+01  6.013e-03 -5.528e+00  6.708e+00  2.518e+01 -4.039e+00]
  DE -ln(L) -6018.5
Starting local hyperparameter optimisation
  Local minimum found at: [-6.000e+00  9.701e-01  1.247e+01  5.981e-03 -5.584e+00  6.384e+00
  2.633e+01 -4.039e+00]
```

]
#block[
```
/Users/benpope/opt/anaconda3/envs/lk/lib/python3.10/site-packages/k2sc/gp.py:98: DeprecationWarning: `np.bool` is a deprecated alias for the builtin `bool`. To silence this warning, use `bool` by itself. Doing this will not modify any behavior and is safe. If you specifically wanted the numpy scalar type, use `np.bool_` here.
Deprecated in NumPy 1.20; for more details and guidance: https://numpy.org/devdocs/release/1.20.0-notes.html#deprecations
  mask = zeros((t1.size,t2.size), np.bool)
```

]
#block[
```
Starting final outlier detection
      8 too high
     31 too low
      0 not finite
Computing time and position trends
  CDPP - raw - 59.783
  CDPP - position component removed - 25.251
  CDPP - full reduction - 25.247
Detrending time 303.0227539539337
```

]
]
#block[
```python
fig, (ax1, ax2) = plt.subplots(1,2, figsize=(12.0,4.0))

detrended = lc.corr_flux-lc.tr_time + np.nanmedian(lc.tr_time)
ax1.plot(lc.time.value,lc.flux.value,'.',label="Uncorrected")
ax1.plot(lc.time.value,detrended.value,'.',label="K2SC")
ax1.set_xlabel('BJD')
ax1.set_ylabel('Flux')

im = ax2.scatter(x,y,c=lc.tr_position-1,cmap='coolwarm',s=1)
ax2.set_xlabel('X')
ax2.set_ylabel('Y')
ax2.set_title('Systematics Model')
plt.colorbar(im,label='GP Correction')

# plt.xlabel('BJD')
# plt.ylabel('Flux')
# plt.title('K2-110b',y=1.01)
```

]
#block[
```python
# search the detrended light curve for transits

bls = BoxLeastSquares(lc.time, detrended, dy=lc.flux_err)

periods = np.linspace(5, 20, 1000)*u.day
durations = np.linspace(0.05, 0.2, 10)*u.day
periodogram = bls.power(periods, durations)

best_period = periods[np.argmax(periodogram.depth_snr.value)]
best_duration = periodogram.duration[np.argmax(periodogram.depth_snr.value)]
best_depth = periodogram.depth[np.argmax(periodogram.depth_snr.value)]
t0 = periodogram.transit_time[np.argmax(periodogram.depth_snr.value)]
```

]
```python
plt.plot(periodogram.period.value,periodogram.depth_snr.value)
plt.xlabel('Period [days]') 
plt.ylabel('Signal-to-noise')
plt.xlim(periods.min().value, periods.max().value)
plt.axvline(best_period.value, color='red', linestyle='dashed')
plt.text(best_period.value+0.1, np.max(periodogram.depth_snr.value), f'Best period: {best_period.value:.4g} d')
```

```
Text(13.958858858858859, 125.15465250327411, 'Best period: 13.86 d')
```

#block[
```python
# what does the model look like?
model = bls.model(lc.time, best_period, best_duration, t0)
plt.plot(lc.time.value, detrended.value)
plt.plot(lc.time.value, model, color='red', lw=1)
```

]
```python
# folded
folded_lc = lc.fold(period=best_period, t0=t0)

folded_times = lc.time.value % best_period.value
folded_args = np.argsort(folded_times)

plt.plot(folded_times[folded_args],detrended.value[folded_args]-1,'.k')
plt.plot(folded_times[folded_args],model[folded_args]-1,color='C3',lw=2)
plt.xlim(2,4)
```

#block[
```
/Users/benpope/opt/anaconda3/envs/lk/lib/python3.10/site-packages/astropy/utils/decorators.py:546: LightkurveDeprecationWarning: "t0" was deprecated in version 2.0 and will be removed in a future version. Use argument "epoch_time" instead.
  return function(*args, **kwargs)
```

]
```
(2.0, 4.0)
```

How did we do this in the paper?

#figure([
#box(image("epic212521166.jpeg"))
], caption: figure.caption(
position: bottom, 
[
EPIC 212521166
]), 
kind: "quarto-float-fig", 
supplement: "Figure", 
)


#box(image("epic212521166.jpg"))

Now let’s look at one last dataset - some asteroseismology data from TESS. Let’s look at the Northern Star - Polaris. "As constant as the Northern Star" - but is it really?

```python
starname = 'Polaris'

search = lightkurve.search_lightcurvefile(starname,mission='TESS',exptime=120)
search
```

#block[
```
/var/folders/vx/lm_q_1ld7c13_fbqfscs9n4w0000gq/T/ipykernel_99792/1715134361.py:3: LightkurveDeprecationWarning: The search_lightcurvefile function is deprecated and may be removed in a future version.
        Use search_lightcurve() instead.
  search = lightkurve.search_lightcurvefile(starname,mission='TESS',exptime=120)
```

]
#block[
```python
lc = search[0].download().PDCSAP_FLUX.normalize().remove_nans()
```

#block[
```
/Users/benpope/opt/anaconda3/envs/lk/lib/python3.10/site-packages/lightkurve/io/tess.py:33: ResourceWarning: unclosed file <_io.BufferedReader name='/Users/benpope/.lightkurve-cache/mastDownload/TESS/tess2019331140908-s0019-0000000303256075-0164-s/tess2019331140908-s0019-0000000303256075-0164-s_lc.fits'>
  lc = read_generic_lightcurve(filename, flux_column=flux_column, time_format="btjd")
/var/folders/vx/lm_q_1ld7c13_fbqfscs9n4w0000gq/T/ipykernel_99792/54457171.py:1: LightkurveDeprecationWarning: The PDCSAP_FLUX function is deprecated and may be removed in a future version.
  lc = search[0].download().PDCSAP_FLUX.normalize().remove_nans()
/Users/benpope/opt/anaconda3/envs/lk/lib/python3.10/site-packages/lightkurve/lightcurve.py:1138: LightkurveWarning: The light curve has a negative median flux (-2.23e+07 electron / s); `normalize()` will therefore divide by a negative number and invert the light curve, which is probablynot what you want
  warnings.warn(
```

]
]
```python
lc.plot()
```

```
<AxesSubplot:xlabel='Time - 2457000 [BTJD days]', ylabel='Normalized Flux'>
```

This is a typical light curve for a Cepheid variable! If we want to determine its period, we can use the Lomb-Scargle periodogram.

This is defined as:

$ P (f) = frac(1, 2 sigma^2) ((sum_i y_i cos (2 pi f t_i))^2 + (sum_i y_i sin (2 pi f t_i))^2) $

where $y_i$ are the data points, $t_i$ are the times, and $sigma$ is the standard deviation of the data.

This is an estimator for the mod-squared Fourier transform suitable for unevenly sampled data. In practice, this is implemented as a linear regression problem with a design matrix that includes sine and cosine terms at a single frequency for all times, and a mean term (and indeed polynomial trends or instrumental systematics models).

Design matrix $X$:

$ X = mat(delim: "[", 1, cos (2 pi f t_1), sin (2 pi f t_1); 1, cos (2 pi f t_2), sin (2 pi f t_2); dots.v, dots.v, dots.v; 1, cos (2 pi f t_N), sin (2 pi f t_N)) $

We then solve this with least squares.

#block[
```python
# lomb-scargle
ls = LombScargle(lc.time, lc.flux, dy=lc.flux_err)
periods = np.linspace(1, 10, 1000)*u.day
frequency = 1/periods
power = ls.power(1/periods)
```

]
```python
plt.plot(periods, power)
plt.xlabel('Period [days]')
plt.ylabel('Power')
plt.xlim(periods.min().value, periods.max().value)
plt.axvline(periods[np.argmax(power)].value, color='red', linestyle='dashed')
plt.text(periods[np.argmax(power)].value+0.1, np.max(power), f'Best period: {periods[np.argmax(power)].value:.4g} d')
```

```
Text(4.063963963963964, 0.9593880168051868, 'Best period: 3.964 d')
```
