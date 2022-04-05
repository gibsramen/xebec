// https://stackoverflow.com/a/53660837
function median(numbers) {
    const sorted = numbers.slice().sort((a, b) => a - b);
    const middle = Math.floor(sorted.length / 2);

    if (sorted.length % 2 === 0) {
        return (sorted[middle - 1] + sorted[middle]) / 2;
    }

    return sorted[middle];
}

// https://stackoverflow.com/a/55297611
const asc = arr => arr.sort((a, b) => a - b);
function quantile(arr, q) {
    const sorted = asc(arr);
    const pos = (sorted.length - 1) * q;
    const base = Math.floor(pos);
    const rest = pos - base;
    if (sorted[base + 1] !== undefined) {
        return sorted[base] + rest * (sorted[base + 1] - sorted[base]);
    } else {
        return sorted[base];
    }
}

const data = big_source.data;
const colArray = data['column'];

const colPhyloArray = [];  // phylogenetic
const colDivMetricArray = [];  // diversity_metric
const colESArray = [];  // effect_size
const colCompArray = []; // comparison
const colColorArray = []; // color
const colMarkerArray = []; // marker

const compCompObj = {};

for (let i = 0; i < colArray.length; i++) {
    if (colArray[i] == cb_obj.value) {
        let phylo = data['phylogenetic'][i];
        let comp = data['comparison'][i];
        let divMetric = data['diversity_metric'][i];
        let effectSize = data['effect_size'][i];
        let color = data['color'][i];
        let marker = data['marker'][i];

        if (comp in compCompObj) {
            compCompObj[comp].push(effectSize)
        } else {
            compCompObj[comp] = [effectSize];
        }

        colPhyloArray.push(phylo);
        colCompArray.push(comp);
        colDivMetricArray.push(divMetric);
        colESArray.push(effectSize);
        colColorArray.push(color);
        colMarkerArray.push(marker);
    } else {
    }
}

const uniqueComps = [...new Set(colCompArray)];

const newBoxSource = {
    'comparison': [],
    'index': [],
    'upper': [],
    'lower': [],
    'q1': [],
    'q2': [],
    'q3': [],
    'qmax': [],
    'qmin': [],
};

var i = 1;
for (let [key, values] of Object.entries(compCompObj)) {
    newBoxSource['comparison'].push(key);
    let qmax = Math.max(...values);
    let qmin = Math.min(...values);
    let q1 = quantile(values, 0.25);
    let q3 = quantile(values, 0.75);
    newBoxSource['qmax'].push(qmax);
    newBoxSource['qmin'].push(qmin);
    newBoxSource['q1'].push(q1);
    newBoxSource['q2'].push(median(values));
    newBoxSource['q3'].push(q3);
    newBoxSource['upper'].push(Math.max(q3, qmax));
    newBoxSource['lower'].push(Math.min(q1, qmin));

    newBoxSource['index'].push(i);
    i = i + 1;
}

box_source.data = newBoxSource;
col_source.data = {
    'phylogenetic': colPhyloArray,
    'diversity_metric': colDivMetricArray,
    'color': colColorArray,
    'marker': colMarkerArray,
    'effect_size': colESArray,
    'comparison': colCompArray
};
// Need to update y_range and factors
yr.end = i - 1;
yr.factors = newBoxSource['comparison'];
box_source.change.emit();
col_source.change.emit();
