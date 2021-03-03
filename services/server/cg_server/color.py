warmColors = [
    # autumn palette, [0, 1, 0.15]
    # '#ff0000',
    # '#ff2600',
    # '#ff4c00',
    # '#ff7300',
    # '#ff9900',
    # '#ffc000',
    # '#ffe600',
    # SHUFFLED:
    "#ff0000",
    "#ff7300",
    "#ff2600",
    "#ff9900",
    "#ff4c00",
    "#ffc000",
    "#ffe600",
]

coolColors = [
    # bwr palette, [0, 0.4, 0.07]
    # '#0000ff',
    # '#2222ff',
    # '#4646ff',
    # '#6969ff',
    # '#8e8eff',
    # '#b2b2ff',
    # cool palette, [0.2, 0.8, 0.05]
    # '#33ccff',
    # '#40bfff',
    # '#4cb3ff',
    # '#59a6ff',
    # '#6699ff',
    # '#738cff',
    # '#7f80ff',
    # '#8c72ff',
    # '#9966ff',
    # '#a659ff',
    # '#b34cff',
    # '#bf40ff',
    # '#cc32ff',
    # SHUFFLED:
    "#0000ff",
    "#6969ff",
    "#2222ff",
    "#8e8eff",
    "#4646ff",
    "#b2b2ff",
    "#33ccff",
    "#7f80ff",
    "#40bfff",
    "#8c72ff",
    "#4cb3ff",
    "#9966ff",
    "#59a6ff",
    "#a659ff",
    "#6699ff",
    "#b34cff",
    "#738cff",
    "#bf40ff",
]

# from https://personal.sron.nl/~pault/#sec:qualitative
# 'vibrant', then 'muted'
snv_colors = [
    "#0077bb",
    "#33bbee",
    "#009988",
    "#ee7733",
    "#cc3311",
    "#ee3377",
    "#332288",
    "#88ccee",
    "#44aa99",
    "#117733",
    "#999933",
    "#ddcc77",
    "#cc6677",
    "#882255",
    "#aa4499",
]

# from https://personal.sron.nl/~pault/#sec:qualitative
# 'vibrant' scheme
clade_colors = [
    "#0077bb",
    "#33bbee",
    "#009988",
    "#ee7733",
    "#cc3311",
    "#ee3377",
]

reds = [
    "#FFF5F0",
    "#FEF1EB",
    "#FEEEE6",
    "#FEEAE1",
    "#FEE7DC",
    "#FEE3D7",
    "#FEE0D2",
    "#FDDACB",
    "#FDD4C3",
    "#FDCEBB",
    "#FCC8B3",
    "#FCC2AB",
    "#FCBCA3",
    "#FCB59B",
    "#FCAF93",
    "#FCA88B",
    "#FCA184",
    "#FC9B7C",
    "#FC9474",
    "#FB8D6D",
    "#FB8767",
]


def get_categorical_colormap(groups):
    # Initialize color map, with the "Other" group as grey by default
    colormap = dict()
    colormap["OTHER"] = "#AAA"

    for i, group in enumerate(groups):
        colormap[group] = clade_colors[i % len(clade_colors)]

    return colormap


def get_snv_colors(snvs):

    colors = []
    for i, snv in enumerate(snvs):
        colors.append(snv_colors[i % len(snv_colors)])

    return colors
