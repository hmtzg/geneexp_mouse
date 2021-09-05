# A treemap R script produced by the REVIGO server at http://revigo.irb.hr/
# If you found REVIGO useful in your work, please cite the following reference:
# Supek F et al. "REVIGO summarizes and visualizes long lists of Gene Ontology
# terms" PLoS ONE 2011. doi:10.1371/journal.pone.0021800

# author: Anton Kratz <anton.kratz@gmail.com>, RIKEN Omics Science Center, Functional Genomics Technology Team, Japan
# created: Fri, Nov 02, 2012  7:25:52 PM
# last change: Fri, Nov 09, 2012  3:20:01 PM

# -----------------------------------------------------------------------------
# If you don't have the treemap package installed, uncomment the following line:
# install.packages( "treemap" );
library(treemap) 								# treemap package by Martijn Tennekes

# Set the working directory if necessary
# setwd("C:/Users/username/workingdir");

# --------------------------------------------------------------------------
# Here is your data from REVIGO. Scroll down for plot configuration options.

revigo.names <- c("term_ID","description","freqInDbPercent","value","uniqueness","dispensability","representative");
revigo.data <- rbind(c("GO:0006732","(obsolete) coenzyme metabolic process",0.395,1.446,1.000,0.000,"(obsolete) coenzyme metabolic process"),
c("GO:0006733","(obsolete) oxidoreduction coenzyme metabolic process",0.395,1.724,1.000,0.000,"(obsolete) oxidoreduction coenzyme metabolic process"),
c("GO:0009108","(obsolete) coenzyme biosynthetic process",0.395,1.587,1.000,0.000,"(obsolete) coenzyme biosynthetic process"),
c("GO:0010257","NADH dehydrogenase complex assembly",0.277,2.195,0.958,0.000,"NADH dehydrogenase complex assembly"),
c("GO:0030198","extracellular matrix organization",1.284,1.808,0.966,0.198,"NADH dehydrogenase complex assembly"),
c("GO:0043062","extracellular structure organization",1.293,1.818,0.966,0.234,"NADH dehydrogenase complex assembly"),
c("GO:0007005","mitochondrion organization",2.085,1.783,0.961,0.247,"NADH dehydrogenase complex assembly"),
c("GO:0051262","protein tetramerization",0.410,1.793,0.956,0.492,"NADH dehydrogenase complex assembly"),
c("GO:0033108","mitochondrial respiratory chain complex assembly",0.449,2.158,0.946,0.574,"NADH dehydrogenase complex assembly"),
c("GO:0032981","mitochondrial respiratory chain complex I assembly",0.277,2.195,0.948,0.671,"NADH dehydrogenase complex assembly"),
c("GO:0035162","embryonic hemopoiesis",0.134,1.845,0.938,0.000,"embryonic hemopoiesis"),
c("GO:0098868","bone growth",0.148,1.787,0.977,0.139,"embryonic hemopoiesis"),
c("GO:0072359","circulatory system development",4.604,1.797,0.963,0.252,"embryonic hemopoiesis"),
c("GO:0048514","blood vessel morphogenesis",2.142,1.789,0.937,0.354,"embryonic hemopoiesis"),
c("GO:0050818","regulation of coagulation",0.339,1.953,0.937,0.000,"regulation of coagulation"),
c("GO:0030195","negative regulation of blood coagulation",0.219,2.077,0.841,0.293,"regulation of coagulation"),
c("GO:1901342","regulation of vasculature development",1.422,1.932,0.912,0.348,"regulation of coagulation"),
c("GO:0046627","negative regulation of insulin receptor signaling pathway",0.191,1.808,0.946,0.381,"regulation of coagulation"),
c("GO:0008360","regulation of cell shape",0.740,1.736,0.917,0.399,"regulation of coagulation"),
c("GO:0045596","negative regulation of cell differentiation",3.512,1.340,0.920,0.481,"regulation of coagulation"),
c("GO:1903034","regulation of response to wounding",0.811,1.629,0.946,0.481,"regulation of coagulation"),
c("GO:1904018","positive regulation of vasculature development",0.802,1.942,0.887,0.545,"regulation of coagulation"),
c("GO:1900046","regulation of hemostasis",0.329,1.863,0.936,0.674,"regulation of coagulation"),
c("GO:0051186","(obsolete) cofactor metabolic process",0.395,1.558,1.000,0.000,"(obsolete) cofactor metabolic process"),
c("GO:0051187","(obsolete) cofactor catabolic process",0.395,1.771,1.000,0.000,"(obsolete) cofactor catabolic process"),
c("GO:0051188","(obsolete) cofactor biosynthetic process",0.395,1.535,1.000,0.000,"(obsolete) cofactor biosynthetic process"),
c("GO:0055114","(obsolete) oxidation-reduction process",0.395,1.788,1.000,0.000,"(obsolete) oxidation-reduction process"),
c("GO:0072577","endothelial cell apoptotic process",0.024,1.748,0.990,0.009,"endothelial cell apoptotic process"),
c("GO:0048659","smooth muscle cell proliferation",0.029,1.585,0.990,0.009,"smooth muscle cell proliferation"),
c("GO:0035924","cellular response to vascular endothelial growth factor stimulus",0.191,1.896,0.916,0.010,"cellular response to vascular endothelial growth factor stimulus"),
c("GO:0048010","vascular endothelial growth factor receptor signaling pathway",0.119,1.818,0.938,0.103,"cellular response to vascular endothelial growth factor stimulus"),
c("GO:0001101","response to acid chemical",0.673,1.469,0.952,0.187,"cellular response to vascular endothelial growth factor stimulus"),
c("GO:0032495","response to muramyl dipeptide",0.095,1.792,0.939,0.247,"cellular response to vascular endothelial growth factor stimulus"),
c("GO:0034097","response to cytokine",4.180,1.544,0.926,0.347,"cellular response to vascular endothelial growth factor stimulus"),
c("GO:0071346","cellular response to interferon-gamma",0.515,1.659,0.872,0.387,"cellular response to vascular endothelial growth factor stimulus"),
c("GO:0000302","response to reactive oxygen species",0.935,1.529,0.928,0.399,"cellular response to vascular endothelial growth factor stimulus"),
c("GO:0007167","enzyme linked receptor protein signaling pathway",2.767,1.441,0.920,0.408,"cellular response to vascular endothelial growth factor stimulus"),
c("GO:0070848","response to growth factor",2.362,1.423,0.931,0.473,"cellular response to vascular endothelial growth factor stimulus"),
c("GO:0071417","cellular response to organonitrogen compound",2.744,1.478,0.891,0.545,"cellular response to vascular endothelial growth factor stimulus"),
c("GO:0034341","response to interferon-gamma",0.630,1.782,0.886,0.589,"cellular response to vascular endothelial growth factor stimulus"),
c("GO:0007178","transmembrane receptor protein serine/threonine kinase signaling pathway",0.921,1.526,0.926,0.611,"cellular response to vascular endothelial growth factor stimulus"),
c("GO:0071560","cellular response to transforming growth factor beta stimulus",0.673,1.549,0.905,0.660,"cellular response to vascular endothelial growth factor stimulus"),
c("GO:0045087","innate immune response",4.213,1.416,0.894,0.672,"cellular response to vascular endothelial growth factor stimulus"),
c("GO:0009060","aerobic respiration",0.649,2.178,0.853,0.011,"aerobic respiration"),
c("GO:0016999","antibiotic metabolic process",0.010,2.070,0.924,0.106,"aerobic respiration"),
c("GO:0019674","NAD metabolic process",0.095,1.807,0.933,0.106,"aerobic respiration"),
c("GO:0006734","NADH metabolic process",0.134,1.863,0.932,0.109,"aerobic respiration"),
c("GO:0032963","collagen metabolic process",0.267,2.019,0.927,0.116,"aerobic respiration"),
c("GO:0140053","mitochondrial gene expression",0.406,2.117,0.912,0.121,"aerobic respiration"),
c("GO:0042744","hydrogen peroxide catabolic process",0.119,1.855,0.914,0.128,"aerobic respiration"),
c("GO:0046034","ATP metabolic process",0.892,1.975,0.919,0.131,"aerobic respiration"),
c("GO:0009161","ribonucleoside monophosphate metabolic process",0.253,1.984,0.710,0.137,"aerobic respiration"),
c("GO:0017144","drug metabolic process",0.339,1.980,0.908,0.141,"aerobic respiration"),
c("GO:0072593","reactive oxygen species metabolic process",0.472,1.520,0.905,0.145,"aerobic respiration"),
c("GO:0006091","generation of precursor metabolites and energy",1.646,1.927,0.893,0.165,"aerobic respiration"),
c("GO:0072524","pyridine-containing compound metabolic process",0.167,1.591,0.849,0.307,"aerobic respiration"),
c("GO:0000959","mitochondrial RNA metabolic process",0.181,1.876,0.844,0.319,"aerobic respiration"),
c("GO:0072350","tricarboxylic acid metabolic process",0.076,1.855,0.817,0.371,"aerobic respiration"),
c("GO:0072521","purine-containing compound metabolic process",1.579,1.768,0.816,0.375,"aerobic respiration"),
c("GO:0046939","nucleotide phosphorylation",0.324,1.646,0.829,0.392,"aerobic respiration"),
c("GO:0032543","mitochondrial translation",0.262,1.970,0.813,0.398,"aerobic respiration"),
c("GO:0006101","citrate metabolic process",0.029,1.871,0.830,0.412,"aerobic respiration"),
c("GO:0072522","purine-containing compound biosynthetic process",0.697,1.703,0.775,0.499,"aerobic respiration"),
c("GO:0043648","dicarboxylic acid metabolic process",0.434,1.777,0.790,0.500,"aerobic respiration"),
c("GO:0055086","nucleobase-containing small molecule metabolic process",2.233,1.706,0.751,0.529,"aerobic respiration"),
c("GO:0046496","nicotinamide nucleotide metabolic process",0.143,1.620,0.732,0.621,"aerobic respiration"),
c("GO:0019693","ribose phosphate metabolic process",1.517,1.743,0.755,0.624,"aerobic respiration"),
c("GO:0046513","ceramide biosynthetic process",0.272,1.682,0.820,0.626,"aerobic respiration"),
c("GO:0006520","cellular amino acid metabolic process",1.169,1.726,0.767,0.638,"aerobic respiration"),
c("GO:0009205","purine ribonucleoside triphosphate metabolic process",0.291,1.972,0.702,0.658,"aerobic respiration"),
c("GO:0009123","nucleoside monophosphate metabolic process",0.348,1.940,0.716,0.676,"aerobic respiration"),
c("GO:1901605","alpha-amino acid metabolic process",0.830,1.712,0.772,0.683,"aerobic respiration"),
c("GO:0046390","ribose phosphate biosynthetic process",0.711,1.665,0.742,0.694,"aerobic respiration"),
c("GO:0006082","organic acid metabolic process",4.013,1.522,0.749,0.696,"aerobic respiration"),
c("GO:0044283","small molecule biosynthetic process",2.243,1.412,0.763,0.697,"aerobic respiration"),
c("GO:1990542","mitochondrial transmembrane transport",0.415,2.122,0.905,0.012,"mitochondrial transmembrane transport"),
c("GO:0043542","endothelial cell migration",0.315,1.841,0.852,0.189,"mitochondrial transmembrane transport"),
c("GO:0070585","protein localization to mitochondrion",0.382,1.906,0.924,0.192,"mitochondrial transmembrane transport"),
c("GO:0006839","mitochondrial transport",0.911,1.717,0.934,0.278,"mitochondrial transmembrane transport"),
c("GO:0006911","phagocytosis, engulfment",1.107,1.699,0.901,0.307,"mitochondrial transmembrane transport"),
c("GO:0015985","energy coupled proton transport, down electrochemical gradient",0.095,1.813,0.915,0.433,"mitochondrial transmembrane transport"),
c("GO:0006816","calcium ion transport",1.121,1.440,0.915,0.534,"mitochondrial transmembrane transport"),
c("GO:0044743","protein transmembrane import into intracellular organelle",0.167,1.846,0.903,0.560,"mitochondrial transmembrane transport"),
c("GO:0050900","leukocyte migration",1.035,1.834,0.858,0.584,"mitochondrial transmembrane transport"),
c("GO:0097529","myeloid leukocyte migration",0.596,1.854,0.846,0.622,"mitochondrial transmembrane transport"),
c("GO:0042330","taxis",2.438,1.586,0.878,0.637,"mitochondrial transmembrane transport"),
c("GO:0001667","ameboidal-type cell migration",0.921,1.558,0.872,0.652,"mitochondrial transmembrane transport"),
c("GO:0060326","cell chemotaxis",0.988,1.792,0.811,0.657,"mitochondrial transmembrane transport"),
c("GO:0031589","cell-substrate adhesion",0.849,1.634,0.970,0.013,"cell-substrate adhesion"),
c("GO:0034446","substrate adhesion-dependent cell spreading",0.267,1.685,0.948,0.571,"cell-substrate adhesion"),
c("GO:0007159","leukocyte cell-cell adhesion",0.282,1.499,0.973,0.574,"cell-substrate adhesion"),
c("GO:2000351","regulation of endothelial cell apoptotic process",0.262,1.748,0.970,0.034,"regulation of endothelial cell apoptotic process"),
c("GO:0043066","negative regulation of apoptotic process",4.652,1.313,0.942,0.537,"regulation of endothelial cell apoptotic process"),
c("GO:0090066","regulation of anatomical structure size",2.763,1.451,0.950,0.035,"regulation of anatomical structure size"),
c("GO:0055065","metal ion homeostasis",3.211,1.437,0.936,0.383,"regulation of anatomical structure size"),
c("GO:0002685","regulation of leukocyte migration",1.059,1.725,0.891,0.043,"regulation of leukocyte migration"),
c("GO:0051924","regulation of calcium ion transport",1.365,1.430,0.928,0.384,"regulation of leukocyte migration"),
c("GO:0002690","positive regulation of leukocyte chemotaxis",0.477,1.706,0.860,0.596,"regulation of leukocyte migration"),
c("GO:0010594","regulation of endothelial cell migration",0.825,1.653,0.858,0.631,"regulation of leukocyte migration"),
c("GO:0002684","positive regulation of immune system process",4.876,1.366,0.935,0.645,"regulation of leukocyte migration"),
c("GO:0051272","positive regulation of cellular component movement",2.939,1.638,0.884,0.694,"regulation of leukocyte migration"),
c("GO:0031343","positive regulation of cell killing",0.324,1.704,0.962,0.044,"positive regulation of cell killing"),
c("GO:0045785","positive regulation of cell adhesion",2.114,1.660,0.943,0.190,"positive regulation of cell killing"),
c("GO:0048660","regulation of smooth muscle cell proliferation",0.825,1.648,0.976,0.048,"regulation of smooth muscle cell proliferation"),
c("GO:0090287","regulation of cellular response to growth factor stimulus",1.412,1.612,0.953,0.051,"regulation of cellular response to growth factor stimulus"),
c("GO:0032101","regulation of response to external stimulus",4.385,1.350,0.949,0.367,"regulation of cellular response to growth factor stimulus"),
c("GO:0030155","regulation of cell adhesion",3.555,1.549,0.972,0.059,"regulation of cell adhesion"),
c("GO:0009584","detection of visible light",0.196,1.733,0.957,0.074,"detection of visible light"),
c("GO:0048265","response to pain",0.191,1.793,0.947,0.098,"response to pain"),
c("GO:0050817","coagulation",0.573,1.634,0.981,0.111,"response to pain"),
c("GO:0090130","tissue migration",0.463,1.586,0.981,0.120,"response to pain"),
c("GO:0009611","response to wounding",1.870,1.469,0.954,0.294,"response to pain"),
c("GO:0042060","wound healing",1.336,1.566,0.951,0.363,"response to pain"));

stuff <- data.frame(revigo.data);
names(stuff) <- revigo.names;

stuff$value <- as.numeric( as.character(stuff$value) );
stuff$freqInDbPercent <- as.numeric( as.character(stuff$freqInDbPercent) );
stuff$uniqueness <- as.numeric( as.character(stuff$uniqueness) );
stuff$dispensability <- as.numeric( as.character(stuff$dispensability) );

saveRDS(stuff, file = 'results/figure4/signifGO_DiCo_REVIGO.rds')

# by default, outputs to a PDF file
# pdf( file="revigo_treemap.pdf", width=16, height=9 ) # width and height are in inches
# 
# # check the tmPlot command documentation for all possible parameters - there are a lot more
# treemap(
#   stuff,
#   index = c("representative","description"),
#   vSize = "value",
#   type = "categorical",
#   vColor = "representative",
#   title = "REVIGO TreeMap",
#   inflate.labels = FALSE,      # set this to TRUE for space-filling group labels - good for posters
#   lowerbound.cex.labels = 0,   # try to draw as many labels as possible (still, some small squares may not get a label)
#   bg.labels = "#CCCCCCAA",   # define background color of group labels
# 								 # "#CCCCCC00" is fully transparent, "#CCCCCCAA" is semi-transparent grey, NA is opaque
#   position.legend = "none"
# )
# 
# dev.off()

