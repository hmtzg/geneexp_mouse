# A treemap R script produced by the Revigo server at http://revigo.irb.hr/
# If you found Revigo useful in your work, please cite the following reference:
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
# Here is your data from Revigo. Scroll down for plot configuration options.

revigo.names <- c("term_ID","description","frequency","value","uniqueness","dispensability","representative");
revigo.data <- rbind(c("GO:0010257","NADH dehydrogenase complex assembly",0.274868489645041,2.19454441909736,0.957208818482524,0,"NADH dehydrogenase complex assembly"),
c("GO:0030198","extracellular matrix organization",1.27956021041657,1.80766465385383,0.965479750791071,0.19698999,"NADH dehydrogenase complex assembly"),
c("GO:0043062","extracellular structure organization",1.28903843419743,1.81815058050673,0.965455006617581,0.23194094,"NADH dehydrogenase complex assembly"),
c("GO:0007005","mitochondrion organization",2.06625278422824,1.78343555633583,0.959777504121228,0.24544534,"NADH dehydrogenase complex assembly"),
c("GO:0051262","protein tetramerization",0.417041846357992,1.7929069773273,0.954641951608306,0.48980078,"NADH dehydrogenase complex assembly"),
c("GO:0033108","mitochondrial respiratory chain complex assembly",0.445476517700583,2.15790309546566,0.944903730624373,0.57115618,"NADH dehydrogenase complex assembly"),
c("GO:0032981","mitochondrial respiratory chain complex I assembly",0.274868489645041,2.19454441909736,0.946963312472648,0.67181885,"NADH dehydrogenase complex assembly"),
c("GO:0035162","embryonic hemopoiesis",0.13743424482252,1.84537086318396,0.93545707685212,0,"embryonic hemopoiesis"),
c("GO:0048265","response to pain",0.184825363726838,1.79331463534831,0.944778685576395,0.10089576,"embryonic hemopoiesis"),
c("GO:0050817","coagulation",0.568693426851808,1.6339347981691,0.979608402696735,0.11339602,"embryonic hemopoiesis"),
c("GO:0090130","tissue migration",0.459693853371878,1.58611288948557,0.980007883815054,0.12318764,"embryonic hemopoiesis"),
c("GO:0098868","bone growth",0.146912468603384,1.78701006957076,0.976170017233147,0.14223856,"embryonic hemopoiesis"),
c("GO:0072359","circulatory system development",4.60167764560921,1.79708920290294,0.961116882596366,0.25276333,"embryonic hemopoiesis"),
c("GO:0009611","response to wounding",1.87194919672053,1.46949189964295,0.952253277194184,0.29298498,"embryonic hemopoiesis"),
c("GO:0048514","blood vessel morphogenesis",2.12786123880385,1.78896528879938,0.93445211469823,0.35270882,"embryonic hemopoiesis"),
c("GO:0042060","wound healing",1.33642955310175,1.56567523673321,0.94925861064201,0.3628753,"embryonic hemopoiesis"),
c("GO:0050818","regulation of coagulation",0.336476944220653,1.95348802335826,0.935654375738262,0,"regulation of coagulation"),
c("GO:0030195","negative regulation of blood coagulation",0.21799914695986,2.07745555722832,0.835958302661788,0.29124771,"regulation of coagulation"),
c("GO:1901342","regulation of vasculature development",1.42647267901995,1.93187631469603,0.909192936425811,0.34626379,"regulation of coagulation"),
c("GO:0046627","negative regulation of insulin receptor signaling pathway",0.189564475617269,1.80785382109806,0.944125736770452,0.38077814,"regulation of coagulation"),
c("GO:0008360","regulation of cell shape",0.748779678688214,1.73602277485749,0.914134833640048,0.39876195,"regulation of coagulation"),
c("GO:1903034","regulation of response to wounding",0.810388133263826,1.62892398621136,0.9446177622291,0.47623125,"regulation of coagulation"),
c("GO:0045596","negative regulation of cell differentiation",3.50220368702905,1.33955966299246,0.917670803036388,0.47968697,"regulation of coagulation"),
c("GO:1904018","positive regulation of vasculature development",0.810388133263826,1.94165360157388,0.883171072583393,0.5433051,"regulation of coagulation"),
c("GO:1900046","regulation of hemostasis",0.32699872043979,1.86267573298743,0.934133513660326,0.67380759,"regulation of coagulation"),
c("GO:0072577","endothelial cell apoptotic process",0.0236955594521587,1.74807445470016,0.989889473531299,0.0085593,"endothelial cell apoptotic process"),
c("GO:0048659","smooth muscle cell proliferation",0.0331737832330221,1.58515997262291,0.989633529741645,0.00876638,"smooth muscle cell proliferation"),
c("GO:0035924","cellular response to vascular endothelial growth factor stimulus",0.189564475617269,1.89569730204042,0.913779308196803,0.01002242,"cellular response to vascular endothelial growth factor stimulus"),
c("GO:0048010","vascular endothelial growth factor receptor signaling pathway",0.118477797260793,1.81819457572826,0.936401802287487,0.10279673,"cellular response to vascular endothelial growth factor stimulus"),
c("GO:0001101","response to acid chemical",0.672953888441306,1.46870838463165,0.9502292464573,0.18571422,"cellular response to vascular endothelial growth factor stimulus"),
c("GO:0032495","response to muramyl dipeptide",0.0947822378086347,1.79204175079869,0.936817618865899,0.24613512,"cellular response to vascular endothelial growth factor stimulus"),
c("GO:0034097","response to cytokine",4.23676603004597,1.54400961934245,0.923935899176712,0.34532387,"cellular response to vascular endothelial growth factor stimulus"),
c("GO:0071346","cellular response to interferon-gamma",0.507084972276195,1.65876771394802,0.868453405496568,0.38555554,"cellular response to vascular endothelial growth factor stimulus"),
c("GO:0000302","response to reactive oxygen species",0.924126818634188,1.52864571982835,0.926439068355943,0.39692201,"cellular response to vascular endothelial growth factor stimulus"),
c("GO:0007167","enzyme linked receptor protein signaling pathway",2.77238045590256,1.44137321926164,0.916889311170801,0.40588282,"cellular response to vascular endothelial growth factor stimulus"),
c("GO:0070848","response to growth factor",2.35059949765414,1.42325428111073,0.928937565189718,0.47111049,"cellular response to vascular endothelial growth factor stimulus"),
c("GO:0071417","cellular response to organonitrogen compound",2.73446756077911,1.4777512506409,0.888123238988773,0.54276895,"cellular response to vascular endothelial growth factor stimulus"),
c("GO:0034341","response to interferon-gamma",0.63030188142742,1.78167405878433,0.882276948091049,0.58495563,"cellular response to vascular endothelial growth factor stimulus"),
c("GO:0007178","transmembrane receptor protein serine/threonine kinase signaling pathway",0.909909482962893,1.52584640057377,0.923852585116075,0.60956953,"cellular response to vascular endothelial growth factor stimulus"),
c("GO:0071560","cellular response to transforming growth factor beta stimulus",0.672953888441306,1.54925524250852,0.901771901498021,0.66092628,"cellular response to vascular endothelial growth factor stimulus"),
c("GO:0045087","innate immune response",4.2557224776077,1.41580659842889,0.891111800150205,0.67115384,"cellular response to vascular endothelial growth factor stimulus"),
c("GO:0009060","aerobic respiration",0.663475664660443,2.17788013015156,0.848353763204318,0.01117305,"aerobic respiration"),
c("GO:0016999","antibiotic metabolic process",0.00947822378086347,2.0696551883112,0.922339098294285,0.10515546,"aerobic respiration"),
c("GO:0019674","NAD metabolic process",0.0947822378086347,1.80737015950375,0.931380028138056,0.10589428,"aerobic respiration"),
c("GO:0006734","NADH metabolic process",0.132695132932089,1.86324179627382,0.92949593570281,0.10895524,"aerobic respiration"),
c("GO:0032963","collagen metabolic process",0.317520496658926,2.01920998461787,0.924091336318179,0.1177836,"aerobic respiration"),
c("GO:0140053","mitochondrial gene expression",0.417041846357992,2.11686600793401,0.909653036110702,0.12084338,"aerobic respiration"),
c("GO:0042744","hydrogen peroxide catabolic process",0.123216909151225,1.85535907963343,0.911052696409121,0.12817996,"aerobic respiration"),
c("GO:0046034","ATP metabolic process",0.905170371072461,1.97501749077856,0.916399017632418,0.13047751,"aerobic respiration"),
c("GO:0009161","ribonucleoside monophosphate metabolic process",0.279607601535472,1.98361051606465,0.700736502476454,0.13782073,"aerobic respiration"),
c("GO:0017144","drug metabolic process",0.331737832330221,1.98030944699768,0.905778052005695,0.14001781,"aerobic respiration"),
c("GO:0072593","reactive oxygen species metabolic process",0.478650300933605,1.51973287423487,0.902639036756424,0.14497409,"aerobic respiration"),
c("GO:0006091","generation of precursor metabolites and energy",1.61129804274679,1.92712424806971,0.890576424702277,0.16421967,"aerobic respiration"),
c("GO:0072524","pyridine-containing compound metabolic process",0.180086251836406,1.59149291651733,0.843181429389163,0.31089411,"aerobic respiration"),
c("GO:0000959","mitochondrial RNA metabolic process",0.194303587507701,1.87614513745872,0.837901825588929,0.32289428,"aerobic respiration"),
c("GO:0072350","tricarboxylic acid metabolic process",0.0758257902469077,1.8545466663298,0.813011334023764,0.37026033,"aerobic respiration"),
c("GO:0072521","purine-containing compound metabolic process",1.63973271408938,1.7680302910806,0.809124497089703,0.37963358,"aerobic respiration"),
c("GO:0046939","nucleotide phosphorylation",0.322259608549358,1.64631977292942,0.823863936965545,0.39352197,"aerobic respiration"),
c("GO:0032543","mitochondrial translation",0.260651153973745,1.96985062019577,0.806444589716373,0.39442857,"aerobic respiration"),
c("GO:0006101","citrate metabolic process",0.0284346713425904,1.87063055641339,0.82588239525929,0.40569017,"aerobic respiration"),
c("GO:0043648","dicarboxylic acid metabolic process",0.426520070138856,1.7770289763791,0.785068297373268,0.49259499,"aerobic respiration"),
c("GO:0072522","purine-containing compound biosynthetic process",0.706127671674328,1.70346984177366,0.766868895530281,0.49958879,"aerobic respiration"),
c("GO:0055086","nucleobase-containing small molecule metabolic process",2.27951281929766,1.70594894404085,0.744186956760698,0.52358107,"aerobic respiration"),
c("GO:0046496","nicotinamide nucleotide metabolic process",0.142173356712952,1.62028673508749,0.725577632377417,0.61821966,"aerobic respiration"),
c("GO:0046513","ceramide biosynthetic process",0.270129377754609,1.6822036532186,0.813633614568641,0.62478048,"aerobic respiration"),
c("GO:0044283","small molecule biosynthetic process",1.96673143452917,1.41238480329416,0.760069681600835,0.62737832,"aerobic respiration"),
c("GO:0019693","ribose phosphate metabolic process",1.60181981896593,1.74295947402024,0.74683583249442,0.62882588,"aerobic respiration"),
c("GO:0006520","cellular amino acid metabolic process",1.1942561963888,1.7260387790229,0.760540448932233,0.62997665,"aerobic respiration"),
c("GO:0009205","purine ribonucleoside triphosphate metabolic process",0.312781384768494,1.97239751168252,0.693288310555493,0.66007091,"aerobic respiration"),
c("GO:1901605","alpha-amino acid metabolic process",0.829344580825553,1.71153410728373,0.766622035016447,0.67541293,"aerobic respiration"),
c("GO:0009123","nucleoside monophosphate metabolic process",0.374389839344107,1.94001722075119,0.707115522967785,0.67703971,"aerobic respiration"),
c("GO:0006082","organic acid metabolic process",4.14672290412777,1.52224394179028,0.741856255168352,0.69457905,"aerobic respiration"),
c("GO:0046390","ribose phosphate biosynthetic process",0.725084119236055,1.66499198182715,0.734575913328457,0.69495881,"aerobic respiration"),
c("GO:1990542","mitochondrial transmembrane transport",0.412302734467561,2.12161636512465,0.902113603653811,0.01160423,"mitochondrial transmembrane transport"),
c("GO:0043542","endothelial cell migration",0.312781384768494,1.84094780018512,0.847540668221855,0.18686285,"mitochondrial transmembrane transport"),
c("GO:0070585","protein localization to mitochondrion",0.38386806312497,1.90612128405551,0.922192473692058,0.19032487,"mitochondrial transmembrane transport"),
c("GO:0006839","mitochondrial transport",0.919387706743756,1.71748292590912,0.932734158110815,0.27526425,"mitochondrial transmembrane transport"),
c("GO:0006911","phagocytosis, engulfment",1.09473484668973,1.69857743490706,0.898226684282235,0.30446914,"mitochondrial transmembrane transport"),
c("GO:0015985","energy coupled proton transport, down electrochemical gradient",0.0947822378086347,1.81264045500275,0.912282262918765,0.43149068,"mitochondrial transmembrane transport"),
c("GO:0006816","calcium ion transport",1.13738685370362,1.43990772961457,0.912379805505268,0.53654655,"mitochondrial transmembrane transport"),
c("GO:0044743","protein transmembrane import into intracellular organelle",0.165868916165111,1.8459043489511,0.900822912489576,0.55890718,"mitochondrial transmembrane transport"),
c("GO:0050900","leukocyte migration",1.04734372778541,1.83350179878073,0.853604193246992,0.5842332,"mitochondrial transmembrane transport"),
c("GO:0097529","myeloid leukocyte migration",0.606606321975262,1.85443203955643,0.841681953512678,0.62427702,"mitochondrial transmembrane transport"),
c("GO:0042330","taxis",2.44064262357234,1.5860699036265,0.873749491660535,0.6369083,"mitochondrial transmembrane transport"),
c("GO:0001667","ameboidal-type cell migration",0.914648594853324,1.55844404495551,0.868746960528966,0.65198233,"mitochondrial transmembrane transport"),
c("GO:0060326","cell chemotaxis",0.999952608881096,1.79243630659235,0.805317021705781,0.65832619,"mitochondrial transmembrane transport"),
c("GO:0031589","cell-substrate adhesion",0.857779252168144,1.63414633107991,0.969287716772411,0.012474,"cell-substrate adhesion"),
c("GO:0034446","substrate adhesion-dependent cell spreading",0.265390265864177,1.68546509875381,0.946024832503721,0.56753331,"cell-substrate adhesion"),
c("GO:0007159","leukocyte cell-cell adhesion",0.308042272878063,1.49925900461396,0.971703866709879,0.57555727,"cell-substrate adhesion"),
c("GO:2000351","regulation of endothelial cell apoptotic process",0.255912042083314,1.74807445470016,0.969282981249902,0.03391012,"regulation of endothelial cell apoptotic process"),
c("GO:0043066","negative regulation of apoptotic process",4.72963366665087,1.31259535912445,0.940320956570233,0.53311942,"regulation of endothelial cell apoptotic process"),
c("GO:0090066","regulation of anatomical structure size",2.78185867968343,1.45056909964866,0.948385000668527,0.03452202,"regulation of anatomical structure size"),
c("GO:0055065","metal ion homeostasis",3.24155253305531,1.43677793210216,0.934308557865291,0.38720704,"regulation of anatomical structure size"),
c("GO:0002685","regulation of leukocyte migration",1.06156106345671,1.72534632478034,0.887842057032082,0.04319784,"regulation of leukocyte migration"),
c("GO:0010959","regulation of metal ion transport",2.12312212691342,1.49159913301524,0.922501209785856,0.40344073,"regulation of leukocyte migration"),
c("GO:0002690","positive regulation of leukocyte chemotaxis",0.473911189043173,1.70621141674869,0.855831079296353,0.59626074,"regulation of leukocyte migration"),
c("GO:0002684","positive regulation of immune system process",5.00450215629591,1.36645701634875,0.933849924201334,0.6282122,"regulation of leukocyte migration"),
c("GO:0010594","regulation of endothelial cell migration",0.815127245154258,1.65342350841399,0.85396755335457,0.63082689,"regulation of leukocyte migration"),
c("GO:0051272","positive regulation of cellular component movement",2.94298848395811,1.63838460347249,0.880609277411721,0.69407971,"regulation of leukocyte migration"),
c("GO:0031343","positive regulation of cell killing",0.350694279891948,1.70422844862656,0.960771322123178,0.04453244,"positive regulation of cell killing"),
c("GO:0045785","positive regulation of cell adhesion",2.25581725984551,1.65980452610109,0.941145025801976,0.19020058,"positive regulation of cell killing"),
c("GO:0048660","regulation of smooth muscle cell proliferation",0.867257475949007,1.64816708434346,0.97541023866573,0.04887124,"regulation of smooth muscle cell proliferation"),
c("GO:0090287","regulation of cellular response to growth factor stimulus",1.42173356712952,1.61220432783008,0.951521261314668,0.05161674,"regulation of cellular response to growth factor stimulus"),
c("GO:0032101","regulation of response to external stimulus",4.49267807212928,1.34986958257853,0.946914318094917,0.36851505,"regulation of cellular response to growth factor stimulus"),
c("GO:0030155","regulation of cell adhesion",3.69650727453675,1.54933569994943,0.970947151499815,0.06014526,"regulation of cell adhesion"),
c("GO:0009584","detection of visible light",0.194303587507701,1.73310942579638,0.955292766764035,0.07351391,"detection of visible light"));

stuff <- data.frame(revigo.data);
names(stuff) <- revigo.names;

stuff$value <- as.numeric( as.character(stuff$value) );
stuff$frequency <- as.numeric( as.character(stuff$frequency) );
stuff$uniqueness <- as.numeric( as.character(stuff$uniqueness) );
stuff$dispensability <- as.numeric( as.character(stuff$dispensability) );

# by default, outputs to a PDF file
pdf("./results/figure4/revigo_treemap.pdf", width=16, height=9 ) # width and height are in inches
# check the tmPlot command documentation for all possible parameters - there are a lot more
treemap(
  stuff,
  index = c("representative","description"),
  vSize = "value",
  type = "categorical",
  vColor = "representative",
  title = "Revigo TreeMap",
  inflate.labels = FALSE,      # set this to TRUE for space-filling group labels - good for posters
  lowerbound.cex.labels = 0,   # try to draw as many labels as possible (still, some small squares may not get a label)
  bg.labels = "#CCCCCCAA",   # define background color of group labels
								 # "#CCCCCC00" is fully transparent, "#CCCCCCAA" is semi-transparent grey, NA is opaque
  position.legend = "none"
)

dev.off()

rr = stuff[stuff$description==stuff$representative,]
range(rr$uniqueness) # 0.8483538 0.9898895
range(rr$dispensability) #  0.00000000 0.07351391

repr = stuff[stuff$description==stuff$representative,c('term_ID', 'representative')]
saveRDS(repr, './results/figure4/revigo_representative.rds')
