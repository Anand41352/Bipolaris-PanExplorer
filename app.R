# ==============================================================================
#   PANGENOME DASHBOARD v11 — Bipolaris sorokiniana (19 Strains)
#   Targeted Patch Update: Fixed Plots, Tables, and Colored Interactors
#   LATEST UPDATE: Panel-wise customizations & customizable labels system
#
#   ► SET YOUR PATH:
DATA_DIR <- "."

# ==============================================================================
# ► GLOBAL LABEL CUSTOMIZATION SECTION ◄
# ==============================================================================

LABELS <- list(
  home_title        = "Bipolaris PanExplorer: An Integrated Pangenome Analysis Dashboard",
  home_strains      = "Strains",
  home_total_genes  = "Total genes",
  home_total_ogs    = "Total Orthogroups",
  home_core_ogs     = "Core Orthogroups",
  home_single_copy  = "Single-copy OGs",
  home_map_title    = "Geographical Distribution of Strains",
  home_metadata_title = "Genome Assemblies & Characteristics",
  home_busco_title  = "BUSCO-Based Genome Quality Assessment",

  summ_cat_prop_title   = "Pangenome Composition",
  summ_per_strain_title = "Genome-Wise Pangenome Partition Summary",
  summ_strain_stack_title = "Genome-Wise Distribution of Pangenome Partitions",

  og_filter_title = "Filter",
  og_table_title  = "Orthogroups",

  pa_settings_title    = "Settings",
  pa_box_title_prefix  = "Orthogroup Presence–Absence Distribution",

  ann_coverage_title = "Functional Annotation Coverage Across Pangenome Partitions",
  ann_table_title    = "Comprehensive Orthogroup Annotation Table",

  cazy_strain_title  = "Genome-Wise Distribution of CAZyme",
  cazy_heat_title    = "Distribution of CAZyme Families Across Genomes",
  cazy_bar_title     = "Global Abundance of CAZyme Families",
  cazy_genelist_title= "CAZyme Catalog",

  eff_type_title     = "Genome-Wise Distribution of Effector Classes",
  eff_cat_title      = "Pangenome Partitioning of Effector",
  eff_genelist_title = "Effector Gene Catalog",

  phi_strain_title   = "PHI Phenotype Profiles Across Genomes",
  phi_table_title    = "PHI-Based Virulence Gene Catalog",

  cafe_ec_title      = "Genome-Wise Gene Family Expansions and Contractions",
  cafe_cat_title     = "Orthogroup Expansion and Contraction Summary",
  cafe_heat_title    = "Gene Family Copy Number Variation Across Genomes",
  cafe_table_title   = "CAFE-Inferred Orthogroup Annotation",

  sel_summary_title  = "Selection Pressure Landscape of Core Orthogroups",
  sel_methods_title  = "Positive Selection Signals Detected by Individual Methods",
  sel_strong_title   = "Annotated High-Confidence Positively Selected Orthogroup",
  sel_all_title      = "Core Orthogroup Selection Pressure Catalog",

  ppi_control_title  = "Network Analysis Controls",
  ppi_network_title  = "Protein–Protein Interaction Network",
  ppi_table_title    = "PPI Network Node Annotation Catalog",

  overview_venn_title   = "Shared Orthogroups Across Functional Categories",
  overview_func_title   = "Distribution of Functional Orthogroups Across Pangenome Partitions",
  overview_master_title = "Integrated Orthogroup Knowledgebase — all 16,980 OGs",

  dl_title = "Export data"
)

# ==============================================================================
library(shiny);           library(shinydashboard); library(DT)
library(ggplot2);         library(dplyr);           library(data.table)
library(plotly);          library(leaflet);          library(RColorBrewer)
library(shinycssloaders); library(shinyWidgets);    library(tidyr)
library(igraph)
library(visNetwork)

ANN_DIR  <- file.path(DATA_DIR, "annotation")
VIR_DIR  <- file.path(DATA_DIR, "virulence_category")
SEQ_DIR  <- file.path(DATA_DIR, "orthogroups")
CAFE_DIR <- file.path(DATA_DIR, "cafe")
SEL_DIR  <- file.path(DATA_DIR, "selection_pressure")
PPI_DIR  <- file.path(DATA_DIR, "ppi")
PPI_CORE <- file.path(PPI_DIR, "core_interactome")
SUMM_DIR <- file.path(DATA_DIR, "summary")

cat_pal <- c(Core="#1565C0", Softcore="#2E7D32", Shell="#E65100", Cloud="#6A1B9A")

ann_type_pal <- c(
  "CAFE"                                              = "#FF6F00",
  "CAFE;CAZymes"                                      = "#F57F17",
  "CAFE;Effector"                                     = "#E65100",
  "CAFE;PHI"                                          = "#BF360C",
  "CAFE;PHI;Positive_Selection"                       = "#4E342E",
  "CAZymes"                                           = "#1B5E20",
  "CAZymes;Effector"                                  = "#2E7D32",
  "CAZymes;Effector;PHI"                              = "#388E3C",
  "CAZymes;Effector;PHI;Positive_Selection"           = "#43A047",
  "CAZymes;PHI"                                       = "#66BB6A",
  "CAZymes;PHI;Positive_Selection"                    = "#A5D6A7",
  "CAZymes;Positive_Selection"                        = "#00695C",
  "Effector"                                          = "#1565C0",
  "Effector;PHI"                                      = "#1976D2",
  "Effector;PHI;Positive_Selection"                   = "#42A5F5",
  "Effector;Positive_Selection"                       = "#90CAF9",
  "Other proteins"                                    = "#9E9E9E",
  "PHI"                                               = "#6A1B9A",
  "PHI;Positive_Selection"                            = "#AB47BC",
  "Positive_Selection"                                = "#E91E63"
)
community_pal <- c("1"="#1565C0","2"="#2E7D32","3"="#E65100","4"="#6A1B9A")

# ── Strain metadata — expanded with all fields from strain_summary.txt ─────
strain_meta <- data.frame(
  Strain = c("BRIP10943a","BRIP27492a","BS112","BS52","D2","SI","LK93",
             "Shoemaker","Yt-6","ND90Pr","ND93-1","WAI2411","WAI2431",
             "WAI2432","WAI3285","WAI3295","WAI3382","WAI3384","WAI3398"),
  Country = c("Australia","Australia","India","India","India","India",
              "China","China","China","USA","USA",rep("Australia",8)),
  Isolation_Source = c("Mycelia and spores","Mycelia and spores","Leaves","Leaves",
                       "Leaves","Leaves","Mycelia","Leaves","leaves","Leaves","Kernels",
                       "N.A.","N.A.","N.A.","N.A.","N.A.","N.A.","N.A.","N.A."),
  Host = c("Wheat","Barley","Wheat","Wheat","Wheat","Wheat","Wheat","Wheat",
            "Avena nuda (naked oat)","Barley","Barley",
            "Wheat","Wheat","Wheat","Wheat","Wheat","Wheat","Wheat","Wheat"),
  Location = c("Biloela, QLD, Australia","Esperance, Western Australia, Australia",
               "Banaras, Varanasi, India","Shillongani, Assam, India",
               "Pusa, Bihar, India","Dharwad, Karnataka, India",
               "Lankao, China","Gansu, China","Zhangjiakou, Hebei, China",
               "North Dakota, United States","North Dakota, Walsh County, USA",
               "Nyngan, NSW, Australia","Gilgandra, NSW, Australia",
               "North Star, NSW, Australia","Inglestone, QLD, Australia",
               "Moonie, QLD, Australia","Pallamallawa, NSW, Australia",
               "Mungindi, QLD, Australia","Walgett, NSW, Australia"),
  Assembly_Level = c("Chromosome","Chromosome","Contig","contig","contig","contig",
                     "Chromosome","Scaffold","contig","Scaffold","Chromosome",
                     "Chromosome","Chromosome","Chromosome","Chromosome",
                     "Chromosome","Chromosome","Chromosome","Chromosome"),
  Pub_Year = c(2019,2019,2019,2022,2022,2022,2022,2020,2024,2013,2023,
               2025,2025,2025,2025,2025,2025,2025,2025),
  Sequencing_Technology = c("PacBio RSII","Oxford Nanopore MinION",
                             "Illumina, Oxford Nanopore, Ion Torrent",
                             "Illumina HiSeq 1000 platform","Illumina HiSeq 1000 platform",
                             "Illumina HiSeq 1000 platform","Oxford Nanopore PromethION",
                             "Illumina HiSeq","Oxford Nanopore",
                             "Sanger, 454, Illumina","PacBio",
                             "Oxford Nanopore MinION","Oxford Nanopore MinION",
                             "Oxford Nanopore MinION","Oxford Nanopore MinION",
                             "Oxford Nanopore MinION","Oxford Nanopore MinION",
                             "Oxford Nanopore MinION","Oxford Nanopore MinION"),
  NCBI_Accession = c("GCA_008452735.1","GCA_008452725.1","GCA_004329375.1",
                     "PRJNA785882","PRJNA785882","PRJNA785882","GCA_024722535.1",
                     "GCA_013416765.1","GCA_037974985.1","GCA_000338995.1",
                     "GCA_033006315.1","PRJNA1017791","PRJNA1017791","PRJNA1017791",
                     "PRJNA1017791","PRJNA1017791","PRJNA1017791","PRJNA1017791","PRJNA1017791"),
  WGS_Accession = c("SRZH01","SRZG01","RCTM01","N.A.","N.A.","N.A.","N.A.",
                    "WNKQ01","JBBEUP01","AEIN01","JASPLA01",
                    "N.A.","N.A.","N.A.","N.A.","N.A.","N.A.","N.A.","N.A."),
  Size_Mb = c(36.92,35.24,37.38,34.41,39.07,33.91,36.23,34.33,38.2,
              34.41,36.48,36.33,36.34,36.77,36.2,36.59,36.68,36.77,36.6),
  Lat = c(-24.40,-33.86,25.32,26.20,25.98,15.46,34.82,36.06,40.82,
          47.55,48.35,-31.56,-31.71,-28.93,-25.03,-27.72,-29.47,-28.98,-30.03),
  Lon = c(150.51,121.89,82.97,92.93,85.67,75.01,114.82,103.83,114.88,
          -101.00,-97.72,147.19,148.66,151.20,148.57,150.37,150.88,149.07,148.12),
  stringsAsFactors = FALSE)

cpal <- c(Australia="#2196F3",India="#FF9800",China="#F44336",USA="#4CAF50")
strain_meta$color <- cpal[strain_meta$Country]

# ── Safe reader ───────────────────────────────────────────────────────────────
rd <- function(path, ...) {
  if (!file.exists(path)) { message("MISSING: ", path); return(data.table()) }
  tryCatch(fread(path, ..., data.table=TRUE),
    error=function(e){ message("ERR: ",path,": ",e$message); data.table() })
}

message("=== Loading data ===")

cats      <- rd(file.path(DATA_DIR,"gene_category_summary.csv"))
if(nrow(cats)>0){ setnames(cats,1,"Orthogroup"); setnames(cats,2,"Category") }
summ      <- rd(file.path(DATA_DIR,"pangenome_summary.txt"), sep="\t")
pa_pre    <- rd(file.path(DATA_DIR,"presence_absence_matrix.csv"))
busco     <- rd(file.path(DATA_DIR,"busco_summary.csv"))

if(nrow(summ)>0){
  for(old in names(summ)){ lw <- tolower(trimws(old))
    if(lw=="genome")   setnames(summ,old,"Genome")
    if(lw=="core")     setnames(summ,old,"Core")
    if(lw=="shell")    setnames(summ,old,"Shell")
    if(lw=="cloud")    setnames(summ,old,"Cloud")
    if(lw=="softcore") setnames(summ,old,"Softcore")
  }
}

basics_lines <- if(file.exists(file.path(DATA_DIR,"basics_stats.txt")))
  readLines(file.path(DATA_DIR,"basics_stats.txt"),warn=FALSE) else character(0)
basics_df <- if(length(basics_lines)>0){
  parts <- strsplit(basics_lines[nchar(trimws(basics_lines))>0],"\t")
  data.frame(Stat=sapply(parts,`[[`,1),
             Val=sapply(parts,function(x) if(length(x)>1) x[2] else ""),
             stringsAsFactors=FALSE)
} else data.frame(Stat=character(),Val=character())

egg      <- rd(file.path(ANN_DIR,"eggnog_annotations.csv"))
nrdb     <- rd(file.path(ANN_DIR,"nr_annotations.tsv"), sep="\t")
ipro     <- rd(file.path(ANN_DIR,"interproscan_annotations.tsv"),    sep="\t")
swsp     <- rd(file.path(ANN_DIR,"swissprot_annotations.csv"))
ann_summ <- rd(file.path(ANN_DIR,"annotation_summary.txt"),         sep="\t")

cazy         <- rd(file.path(VIR_DIR,"cazyme_annotations.csv"))
cazy_sub_raw <- rd(file.path(VIR_DIR,"cazyme_subclass_distribution.csv"), sep="\t")
cazy_og      <- rd(file.path(VIR_DIR,"cazyme_orthogroups.txt"), sep="\t")
eff          <- rd(file.path(VIR_DIR,"effector_annotations.csv"))
eff_summ     <- rd(file.path(VIR_DIR,"effector_summary.txt"), sep="\t")
phi_strain   <- rd(file.path(VIR_DIR,"phi_phenotype_distribution.csv"))
phi_db       <- rd(file.path(VIR_DIR,"phi_annotations.tsv"), sep="\t")

cafe_ann  <- rd(file.path(CAFE_DIR,"cafe_expanded_contracted_annotations.csv"))
cafe_heat <- rd(file.path(CAFE_DIR,"cafe_heatmap_data.csv"))
cafe_sum  <- rd(file.path(CAFE_DIR,"cafe_strain_summary.csv"))

sel_core   <- rd(file.path(SEL_DIR,"core_positive_selection.csv"))
sel_strong <- rd(file.path(SEL_DIR,"positively_selected_gene_annotations.csv"))

gene_cat_summ  <- rd(file.path(SUMM_DIR,"gene_category_summary.csv"))
cafe_list      <- rd(file.path(SUMM_DIR,"cafe_summary.txt"), header=FALSE)
cazy_list      <- rd(file.path(SUMM_DIR,"cazyme_summary.txt"), header=FALSE)
eff_list       <- rd(file.path(SUMM_DIR,"effector_summary.txt"), header=FALSE)
phi_list       <- rd(file.path(SUMM_DIR,"phi_summary.txt"), header=FALSE)
sel_list       <- rd(file.path(SUMM_DIR,"positive_selection_summary.txt"), header=FALSE)
top10_ppi_list <- rd(file.path(SUMM_DIR,"top10_percent_ppi_nodes.csv"))

top10_nodes  <- rd(file.path(PPI_DIR,"Top10_Nodes.csv"))
ppi_ann_full <- rd(file.path(PPI_CORE,"orthogroup_annotations.csv"))

message("=== Files loaded ===")
fmt <- function(x) format(as.integer(x),big.mark=",",scientific=FALSE)

TOTAL_OG    <- 16980L
TOTAL_GENES <- 257142L
CORE_OG     <- 10323L
SC_OG       <- 8396L

cat_counts <- if(nrow(cats)>0) cats[,.N,by=Category] else
  data.table(Category=c("Core","Softcore","Shell","Cloud"),N=c(10323L,1444L,3350L,1864L))

summ_long <- if(nrow(summ)>0 && all(c("Genome","Core","Shell","Cloud","Softcore") %in% names(summ)))
  melt(summ,id.vars="Genome",measure.vars=c("Core","Shell","Cloud","Softcore"),
       variable.name="Category",value.name="Count") else data.table()

cazy_sub <- if(nrow(cazy_sub_raw)>0){
  d <- copy(cazy_sub_raw); setnames(d,names(d)[1],"Subclass")
  d[,Total:=rowSums(.SD,na.rm=TRUE),.SDcols=setdiff(names(d),"Subclass")]
  d } else data.table()

if(nrow(phi_db)>0){
  if("query" %in% names(phi_db)){
    phi_db[,OG:=sub("\\|.*","",query)]; phi_db[,Gene_ID:=sub(".*\\|","",query)]
  }
  if("sseqid" %in% names(phi_db)){
    phi_db[,PHI_Gene :=sapply(strsplit(as.character(sseqid),"#"),function(x) if(length(x)>=3) x[3] else "")]
    phi_db[,Organism :=sapply(strsplit(as.character(sseqid),"#"),function(x) if(length(x)>=5) x[5] else "")]
    phi_db[,Phenotype:=sapply(strsplit(as.character(sseqid),"#"),function(x) if(length(x)>=6) x[6] else "")]
  }
}

if(nrow(ipro)>0){
  ipro[,OG:=sub("\\|.*","",get(names(ipro)[1]))]
  ipro_best <- ipro[!is.na(InterPro_Accession) &
    nchar(trimws(as.character(InterPro_Accession)))>0 &
    trimws(as.character(InterPro_Accession))!="-",
    .(IPR_Acc=InterPro_Accession[1],IPR_Desc=InterPro_Description[1],GO=GO_Terms[1]),by=OG]
} else ipro_best <- data.table(OG=character(),IPR_Acc=character(),IPR_Desc=character(),GO=character())

if(nrow(swsp)>0 && "query" %in% names(swsp)){
  swsp[,OG:=sub("\\|.*","",query)]
  swsp_best <- swsp[order(-bitscore),.SD[1],by=OG][,.(OG,SP_Hit=target,SP_pident=pident,SP_evalue=evalue,bitscore)]
} else swsp_best <- data.table(OG=character(),SP_Hit=character(),SP_pident=numeric(),SP_evalue=numeric(),bitscore=numeric())

ann_master <- copy(cats)
if(nrow(egg)>0){
  oc <- if("OG" %in% names(egg)) "OG" else if("Orthogroup" %in% names(egg)) "Orthogroup" else names(egg)[1]
  egg_sub <- egg[,.(Orthogroup=get(oc),
    COG=if("COG_category" %in% names(egg)) COG_category else NA_character_,
    EGG_Desc=if("Description" %in% names(egg)) Description else NA_character_,
    EGG_Name=if("Preferred_name" %in% names(egg)) Preferred_name else NA_character_,
    GOs=if("GOs" %in% names(egg)) GOs else NA_character_,
    KEGG=if("KEGG_ko" %in% names(egg)) KEGG_ko else NA_character_,
    PFAMs=if("PFAMs" %in% names(egg)) PFAMs else NA_character_)]
  ann_master <- merge(ann_master,egg_sub,by="Orthogroup",all.x=TRUE)
}
if(nrow(nrdb)>0){
  on <- if("OG" %in% names(nrdb)) "OG" else names(nrdb)[1]
  nr_sub <- nrdb[,.(Orthogroup=get(on),
    NR_Hit=if("NR_Hit" %in% names(nrdb)) NR_Hit else NA_character_,
    NR_Identity=if("NR_Identity" %in% names(nrdb)) NR_Identity else NA_real_,
    NR_Desc=if("NR_Description" %in% names(nrdb)) NR_Description else
            if("NR_Function" %in% names(nrdb)) NR_Function else NA_character_,
    NR_Species=if("NR_Species" %in% names(nrdb)) NR_Species else NA_character_)]
  ann_master <- merge(ann_master,nr_sub,by="Orthogroup",all.x=TRUE)
}
if(nrow(ipro_best)>0) ann_master <- merge(ann_master,ipro_best,by.x="Orthogroup",by.y="OG",all.x=TRUE)
if(nrow(swsp_best)>0) ann_master <- merge(ann_master,swsp_best,by.x="Orthogroup",by.y="OG",all.x=TRUE)

make_ann_vec <- function(list_df, list_col=1){
  if(nrow(list_df)==0) return(character(0))
  as.character(list_df[[list_col]])
}

cafe_ogs  <- make_ann_vec(cafe_list)
cazy_ogs  <- make_ann_vec(cazy_list)
eff_ogs   <- make_ann_vec(eff_list)
phi_ogs   <- make_ann_vec(phi_list)
sel_ogs   <- make_ann_vec(sel_list)
top10_ogs <- if(nrow(top10_ppi_list)>0) as.character(top10_ppi_list[[1]]) else character(0)

addResourcePath(prefix='ppi_mesh',   directoryPath=file.path(DATA_DIR,"ppi"))
addResourcePath(prefix='venn_upset', directoryPath=file.path(DATA_DIR,"summary"))

# =============================================================================
# UI
# =============================================================================
ui <- dashboardPage(
  skin="blue",
  # ── FIX: this `title=` argument sets the actual HTML <head><title> tag
  # (the text shown in the browser tab / bookmark). It was missing before —
  # dashboardHeader(title=...) only controls the navbar text, NOT the HTML
  # document title, so the browser tab was blank/generic ("Shiny").
  title = "Bipolaris PanExplorer | Pangenome Dashboard",
  dashboardHeader(title=tags$span(style="font-size:12px;font-weight:700;color:#fff;",
    "Bipolaris sorokiniana Pangenome"), titleWidth=240),
  dashboardSidebar(width=240,
    sidebarMenu(id="tabs",
      menuItem("Home & World Map",      tabName="home",   icon=icon("earth-americas")),
      menuItem("Pangenome Summary",     tabName="sumtab", icon=icon("chart-pie")),
      menuItem("Orthogroup Explorer",   tabName="og",     icon=icon("magnifying-glass")),
      menuItem("Presence / Absence",    tabName="pa",     icon=icon("table-cells")),
      menuItem("Annotation Hub",        tabName="ann",    icon=icon("tags")),
      menuItem("CAZymes",               tabName="cazy",   icon=icon("atom")),
      menuItem("Effectors",             tabName="eff",    icon=icon("biohazard")),
      menuItem("PHI Database",          tabName="phi",    icon=icon("shield-virus")),
      menuItem("CAFE Analysis",         tabName="cafe",   icon=icon("tree")),
      menuItem("Selection Pressure",    tabName="sel",    icon=icon("arrow-trend-up")),
      menuItem("Core PPI Network",      tabName="ppi",    icon=icon("circle-nodes")),
      menuItem("Pangenome Overview",    tabName="ogsumm", icon=icon("chart-bar")),
      menuItem("Downloads",             tabName="dl",     icon=icon("download"))
    )
  ),
  dashboardBody(
    tags$head(
      # ── FIX: explicit HTML heading block. Some browsers / linters flag a
      # dashboard body that has no semantic <h1> at all; shinydashboard's
      # box/menu titles are all <h3>/<span>, so we add one visually-hidden
      # <h1> for accessibility/SEO while keeping the visual banner (.ms) as-is.
      tags$title("Bipolaris PanExplorer | Pangenome Dashboard"),
      tags$style(HTML("
      body,.content-wrapper{background:#f0f2f5!important}
      .visually-hidden-h1{position:absolute;width:1px;height:1px;padding:0;margin:-1px;
          overflow:hidden;clip:rect(0,0,0,0);white-space:nowrap;border:0}
      .ms{background:linear-gradient(135deg,#0d47a1,#1565c0 55%,#0277bd);color:#fff;
          padding:15px 22px;border-radius:8px;margin-bottom:14px;box-shadow:0 2px 8px rgba(0,0,0,.2)}
      .ms .ttl{font-size:28px;font-weight:700;line-height:1.45;margin:0}
      .ms .subtitle{font-size:14px;font-style:italic;color:#cce4ff;margin-top:6px;line-height:1.4}
      .kpi{background:#fff;border-radius:9px;padding:14px 16px;
           box-shadow:0 1px 4px rgba(0,0,0,.1);border-left:5px solid #1565c0;margin-bottom:3px}
      .kpi.g{border-left-color:#2e7d32}.kpi.o{border-left-color:#e65100}
      .kpi.p{border-left-color:#6a1b9a}.kpi.t{border-left-color:#00695c}
      .kpi .v{font-size:26px;font-weight:700;color:#1a1a2e;line-height:1.1}
      .kpi .l{font-size:11px;color:#666;margin-top:2px}
      .box{border-radius:9px!important}
      .ok{background:#e8f5e9;border-left:4px solid #2e7d32;padding:8px 12px;
          border-radius:5px;font-size:12px;margin-bottom:10px}
      .warn{background:#fff3e0;border-left-color:#e65100}
      .ann-pill{display:inline-block;padding:2px 8px;border-radius:11px;font-size:11px;font-weight:600}
      .pill-core{background:#dbeafe;color:#1565c0}
      .pill-soft{background:#dcfce7;color:#166534}
      .pill-shell{background:#ffedd5;color:#9a3412}
      .pill-cloud{background:#f3e8ff;color:#6b21a8}
      .heat-info{background:#e3f2fd;border-left:4px solid #1565c0;padding:6px 10px;
                 border-radius:5px;font-size:11px;margin-bottom:8px}
      .pager-btn{margin:4px 2px}
      .legend-dot{display:inline-block;width:12px;height:12px;border-radius:50%;
                  margin-right:4px;vertical-align:middle}
    "))),
    # ── FIX: real, semantic <h1> heading block for the whole dashboard.
    # Visually hidden (so your existing banner design is untouched) but
    # present in the DOM — this is the "HTML heading block" that was missing.
    tags$h1(class="visually-hidden-h1", LABELS$home_title),
    tabItems(

      # ══════════════════════════════════════════════════════════════════
      # HOME & WORLD MAP
      # Layout:  ROW1 = KPI cards
      #          ROW2 = World map (full width 12)
      #          ROW3 = Strain metadata table (full width 12)
      #          ROW4 = BUSCO figure (full width 12)
      # Changes:
      #  - Subtitle added (italic, Bipolaris sorokiniana in em)
      #  - BUSCO title removed from box header (kept only in plotly title)
      #  - Strain metadata now shows ALL 10 columns from strain_summary.txt
      #  - World map is full width (no sidebar metadata column)
      # ══════════════════════════════════════════════════════════════════
      tabItem(tabName="home",
        # ── Banner ──────────────────────────────────────────────────────
        tags$div(class="ms",
          tags$p(class="ttl", LABELS$home_title),
          tags$p(class="subtitle",
            "Interactive Exploration of Pangenomic Diversity, Functional Annotation, ",
            "Virulence Factors, and Evolutionary Dynamics, PPI Networks in ",
            tags$em("Bipolaris sorokiniana")
          )
        ),
        uiOutput("status_banner"),

        # ── ROW 1: KPI cards ────────────────────────────────────────────
        fluidRow(
          column(2,tags$div(class="kpi",       tags$div(class="v","19"),          tags$div(class="l",LABELS$home_strains))),
          column(2,tags$div(class="kpi g",     tags$div(class="v",fmt(TOTAL_GENES)),tags$div(class="l",LABELS$home_total_genes))),
          column(2,tags$div(class="kpi o",     tags$div(class="v",fmt(TOTAL_OG)),  tags$div(class="l",LABELS$home_total_ogs))),
          column(2,tags$div(class="kpi p",     tags$div(class="v",fmt(CORE_OG)),   tags$div(class="l",LABELS$home_core_ogs))),
          column(2,tags$div(class="kpi t",     tags$div(class="v",fmt(SC_OG)),     tags$div(class="l",LABELS$home_single_copy)))
        ),br(),

        # ── ROW 2: World map — full width ───────────────────────────────
        fluidRow(
          box(title=tags$span(icon("earth-americas")," ",LABELS$home_map_title),
              width=12, status="primary", solidHeader=TRUE,
              leafletOutput("world_map", height="460px"),
              tags$div(style="margin-top:7px;font-size:12px;color:#555;",
                tags$b(style="color:#2196F3;","● Australia  "),
                tags$b(style="color:#FF9800;","● India  "),
                tags$b(style="color:#F44336;","● China  "),
                tags$b(style="color:#4CAF50;","● USA"),
                "  — Click marker for details"))
        ),

        # ── ROW 3: Strain metadata — full width, all columns ────────────
        fluidRow(
          box(title=tags$span(icon("list")," ",LABELS$home_metadata_title),
              width=12, status="info", solidHeader=TRUE,
              DT::dataTableOutput("tbl_strains"))
        ),

        # ── ROW 4: BUSCO — full width ───────────────────────────────────
        fluidRow(
          box(title=tags$span(icon("chart-bar")," ",LABELS$home_busco_title),
              width=12, status="primary", solidHeader=TRUE,
              plotlyOutput("plt_busco_split", height="650px"))
        )
      ),

      # ══════════════════════════════════════════════════════════════════
      # PANGENOME SUMMARY
      # ══════════════════════════════════════════════════════════════════
      tabItem(tabName="sumtab",
        fluidRow(
          box(title=LABELS$summ_cat_prop_title,width=6,status="primary",solidHeader=TRUE,
              plotlyOutput("plt_cat_donut",height="480px")),
          box(title=LABELS$summ_per_strain_title,width=6,status="info",solidHeader=TRUE,
              DT::dataTableOutput("tbl_summary"))
        ),
        fluidRow(
          box(title=LABELS$summ_strain_stack_title,width=12,status="primary",solidHeader=TRUE,
              plotlyOutput("plt_strain_stack",height="400px"))
        )
      ),

      # ══════════════════════════════════════════════════════════════════
      # ORTHOGROUP EXPLORER
      # ══════════════════════════════════════════════════════════════════
      tabItem(tabName="og",
        fluidRow(
          column(3,box(title=LABELS$og_filter_title,width=NULL,status="primary",solidHeader=TRUE,
            textInput("og_q","Search OG ID","",placeholder="e.g. OG0000042"),
            selectInput("og_cat","Category",c("All","Core","Softcore","Shell","Cloud")),
            actionButton("og_go","Search",icon=icon("magnifying-glass"),
                         class="btn-primary",style="width:100%"),br(),br(),
            tags$small(style="color:#888;","Click a row to view members."))),
          column(9,box(title=LABELS$og_table_title,width=NULL,status="info",solidHeader=TRUE,
            DT::dataTableOutput("tbl_og")))
        ),
        uiOutput("og_detail_ui")
      ),

      # ══════════════════════════════════════════════════════════════════
      # PRESENCE / ABSENCE
      # ══════════════════════════════════════════════════════════════════
      tabItem(tabName="pa",
        fluidRow(
          column(3,box(title=LABELS$pa_settings_title,width=NULL,status="primary",solidHeader=TRUE,
            selectInput("pa_cat","Category",c("Shell","Cloud","Softcore","Core","All"),"Shell"),
            tags$hr(),
            tags$b("Search single OG:"),br(),
            textInput("pa_search","OG ID","",placeholder="e.g. OG0000006"),
            actionButton("pa_search_btn","Search OG",icon=icon("magnifying-glass"),
                         class="btn-info btn-sm",style="width:100%"),br(),br(),
            tags$b("Paginated view:"),
            uiOutput("pa_page_info"),
            fluidRow(
              column(6,actionButton("pa_prev","◀ Prev",class="btn-sm btn-default pager-btn",style="width:100%")),
              column(6,actionButton("pa_next","Next ▶",class="btn-sm btn-primary pager-btn",style="width:100%"))
            )
          )),
          column(9,
            uiOutput("pa_search_result"),
            box(title=uiOutput("pa_box_title"),width=NULL,status="info",solidHeader=TRUE,
              withSpinner(plotlyOutput("plt_pa",height="1200px"))))
        )
      ),

      # ══════════════════════════════════════════════════════════════════
      # ANNOTATION HUB
      # ══════════════════════════════════════════════════════════════════
      tabItem(tabName="ann",
        fluidRow(
          box(title=LABELS$ann_coverage_title,width=12,status="primary",solidHeader=TRUE,
              plotlyOutput("plt_ann_coverage",height="300px"),
              tags$small(style="color:#888;",
                "Legend: Total orthogroups · NRDB · EggNOG · InterProScan · SwissProt"))
        ),
        fluidRow(
          box(title=tags$span(icon("magnifying-glass")," ",LABELS$ann_table_title),
              width=12,status="success",solidHeader=TRUE,
            fluidRow(
              column(3,selectInput("ann_cat_f","Category",c("All","Core","Softcore","Shell","Cloud"))),
              column(3,textInput("ann_q","Search OG / description","",placeholder="gene name, OG ID…")),
              column(2,br(),actionButton("ann_go","Search",icon=icon("magnifying-glass"),class="btn-primary"))
            ),br(),
            DT::dataTableOutput("tbl_ann"))
        )
      ),

      # ══════════════════════════════════════════════════════════════════
      # CAZymes
      # ══════════════════════════════════════════════════════════════════
      tabItem(tabName="cazy",
        fluidRow(
          box(title=LABELS$cazy_strain_title,width=12,status="primary",solidHeader=TRUE,
              plotlyOutput("plt_cazy_strain",height="360px"))
        ),
        fluidRow(
          column(3,box(title="Heatmap controls",width=NULL,status="primary",solidHeader=TRUE,
            textInput("cazy_heat_q","Search family","",placeholder="e.g. GH18"),
            actionButton("cazy_heat_go","Filter / Reset",icon=icon("filter"),
                         class="btn-primary",style="width:100%"),br(),br(),
            tags$b("Pagination (30 families/page):"),
            uiOutput("cazy_heat_page_info"),
            fluidRow(
              column(6,actionButton("cazy_heat_prev","◀",class="btn-sm btn-default",style="width:100%")),
              column(6,actionButton("cazy_heat_next","▶",class="btn-sm btn-primary",style="width:100%"))
            )
          )),
          column(9,
            box(title=LABELS$cazy_heat_title,width=NULL,status="info",solidHeader=TRUE,
              withSpinner(plotlyOutput("plt_cazy_heat",height="680px"))))
        ),
        fluidRow(
          box(title=LABELS$cazy_bar_title,width=12,status="warning",solidHeader=TRUE,
              withSpinner(plotlyOutput("plt_cazy_bar",height="460px")))
        ),
        fluidRow(
          box(title=LABELS$cazy_genelist_title,width=12,status="success",solidHeader=TRUE,
            fluidRow(
              column(3,selectInput("cazy_f_strain","Strain",
                c("All",sort(unique(if(nrow(cazy)>0) cazy$Strain else character()))))),
              column(3,selectInput("cazy_f_cat","Pangenome category",
                c("All","Core","Softcore","Shell","Cloud"))),
              column(3,if(nrow(cazy)>0 && "CAZyme_Subclass" %in% names(cazy))
                selectInput("cazy_f_sub","CAZyme subclass",c("All",sort(unique(cazy$CAZyme_Subclass))))
                else selectInput("cazy_f_sub","CAZyme subclass",c("All")))
            ),br(),
            DT::dataTableOutput("tbl_cazy"))
        )
      ),

      # ══════════════════════════════════════════════════════════════════
      # EFFECTORS
      # ══════════════════════════════════════════════════════════════════
      tabItem(tabName="eff",
        fluidRow(
          box(title=LABELS$eff_type_title,width=8,status="primary",solidHeader=TRUE,
              plotlyOutput("plt_eff_type",height="380px")),
          box(title=LABELS$eff_cat_title,width=4,status="warning",solidHeader=TRUE,
              plotlyOutput("plt_eff_cat",height="380px"))
        ),
        fluidRow(
          box(title=LABELS$eff_genelist_title,width=12,status="success",solidHeader=TRUE,
            fluidRow(
              column(3,selectInput("eff_f_strain","Strain",
                c("All",if(nrow(eff_summ)>0) sort(eff_summ$Strain) else character()))),
              column(3,selectInput("eff_f_cat","Pangenome category",
                c("All","Core","Softcore","Shell","Cloud"))),
              column(3,if(nrow(eff)>0 && "Effector_Type" %in% names(eff))
                selectInput("eff_f_type","Effector type",c("All",sort(unique(eff$Effector_Type))))
                else selectInput("eff_f_type","Effector type",c("All")))
            ),br(),
            DT::dataTableOutput("tbl_eff"))
        )
      ),

      # ══════════════════════════════════════════════════════════════════
      # PHI DATABASE
      # ══════════════════════════════════════════════════════════════════
      tabItem(tabName="phi",
        fluidRow(
          box(title=LABELS$phi_strain_title,width=12,status="primary",solidHeader=TRUE,
              withSpinner(plotlyOutput("plt_phi_strain",height="540px")))
        ),
        fluidRow(
          box(title=LABELS$phi_table_title,width=12,status="success",solidHeader=TRUE,
              DT::dataTableOutput("tbl_phi"))
        )
      ),

      # ══════════════════════════════════════════════════════════════════
      # CAFE
      # ══════════════════════════════════════════════════════════════════
      tabItem(tabName="cafe",
        fluidRow(
          box(title=LABELS$cafe_ec_title,width=7,status="primary",solidHeader=TRUE,
              plotlyOutput("plt_cafe_ec",height="360px")),
          box(title=LABELS$cafe_cat_title,width=5,status="info",solidHeader=TRUE,
              plotlyOutput("plt_cafe_cat",height="360px"))
        ),
        fluidRow(
          column(3,box(title="Heatmap controls",width=NULL,status="primary",solidHeader=TRUE,
            tags$b("Pagination (25 families/page):"),
            uiOutput("cafe_heat_page_info"),
            fluidRow(
              column(6,actionButton("cafe_heat_prev","◀",class="btn-sm btn-default",style="width:100%")),
              column(6,actionButton("cafe_heat_next","▶",class="btn-sm btn-primary",style="width:100%"))
            ),br(),
            tags$small(style="color:#888;","Exact gene-copy values.")
          )),
          column(9,
            box(title=LABELS$cafe_heat_title,width=NULL,status="warning",solidHeader=TRUE,
              withSpinner(plotlyOutput("plt_cafe_heat",height="680px"))))
        ),
        fluidRow(
          box(title=LABELS$cafe_table_title,width=12,status="success",solidHeader=TRUE,
              DT::dataTableOutput("tbl_cafe"))
        )
      ),

      # ══════════════════════════════════════════════════════════════════
      # SELECTION PRESSURE
      # ══════════════════════════════════════════════════════════════════
      tabItem(tabName="sel",
        fluidRow(
          box(title=LABELS$sel_summary_title,width=5,status="primary",solidHeader=TRUE,
              plotlyOutput("plt_sel_summary",height="360px")),
          box(title=LABELS$sel_methods_title,width=7,status="info",solidHeader=TRUE,
              plotlyOutput("plt_sel_methods",height="360px"))
        ),
        fluidRow(
          box(title=LABELS$sel_strong_title,width=12,status="success",solidHeader=TRUE,
              DT::dataTableOutput("tbl_sel_strong"))
        ),
        fluidRow(
          box(title=LABELS$sel_all_title,width=12,status="warning",solidHeader=TRUE,
              DT::dataTableOutput("tbl_sel_all"))
        )
      ),

      # ══════════════════════════════════════════════════════════════════
      # CORE PPI NETWORK
      # ══════════════════════════════════════════════════════════════════
      tabItem(tabName="ppi",
        fluidRow(
          column(3,
            box(title=LABELS$ppi_control_title,width=NULL,status="primary",solidHeader=TRUE,
              tags$div(class="heat-info",
                tags$b("Top 10% Subnetwork:"),
                " Rendering via high-resolution pre-computed topological canvas map.",br(),br(),
                tags$b("Search any OG:"), " Interacting partner nodes color-coded dynamically by multi-omic annotation profiles."),
              selectInput("ppi_color_by","Colour nodes by",
                c("Community (Louvain)"="community","Annotation type"="ann")),
              textInput("ppi_search_og","Search OG ID","",placeholder="e.g. OG0003930"),
              actionButton("ppi_search_btn","Show interactions",
                           icon=icon("circle-nodes"),class="btn-info",style="width:100%"),
              br(),br(),
              actionButton("ppi_reset_btn","Reset to Top 10%",
                           icon=icon("rotate"),class="btn-default",style="width:100%"),
              br(),br(),
              tags$b("Annotation legend:"),br(),
              uiOutput("ppi_legend_ui")
            )
          ),
          column(9,
            box(title=uiOutput("ppi_net_title"),width=NULL,status="info",solidHeader=TRUE,
              withSpinner(uiOutput("ui_ppi_canvas")))
          )
        ),
        fluidRow(
          box(title=LABELS$ppi_table_title,width=12,status="success",solidHeader=TRUE,
            fluidRow(
              column(3,textInput("ppi_tbl_q","Search OG or annotation","",
                                 placeholder="OG ID or description")),
              column(2,br(),actionButton("ppi_tbl_go","Search",class="btn-info btn-sm"))
            ),br(),
            DT::dataTableOutput("tbl_ppi_ann"))
        )
      ),

      # ══════════════════════════════════════════════════════════════════
      # PANGENOME OVERVIEW
      # Changes:
      #  - UpSet plot removed
      #  - Venn diagram: white background, no caption label below image
      #  - Correct SVG filename: venn_diagram.svg
      # ══════════════════════════════════════════════════════════════════
      tabItem(tabName="ogsumm",
        fluidRow(
          box(title=LABELS$overview_venn_title,
              width=12, status="primary", solidHeader=TRUE,
              # White background, image centred, NO caption text below
              tags$div(
                style="display:flex; align-items:center; justify-content:center;
                       background:#ffffff; padding:20px; border-radius:8px;",
                tags$img(
                  src="venn_upset/venn_diagram.svg",
                  style="max-width:100%; max-height:620px; object-fit:contain; display:block;"
                )
              )
          )
        ),
        fluidRow(
          box(title=LABELS$overview_func_title,width=12,status="success",solidHeader=TRUE,
              plotlyOutput("plt_func_cat_dist",height="450px"))
        ),
        fluidRow(
          box(title=LABELS$overview_master_title,width=12,status="primary",solidHeader=TRUE,
            fluidRow(
              column(2,selectInput("master_cat","Category",
                c("All","Core","Softcore","Shell","Cloud"))),
              column(4,textInput("master_q","Search OG / annotation","",
                placeholder="OG ID, gene name...")),
              column(3,selectInput("master_func","Functional",
                c("All","CAFE","CAZymes","Effector","PHI","Selection","PPI","Multi","None"))),
              column(3,br(),actionButton("master_go","Filter",class="btn-info",style="width:100%"))
            ),br(),
            DT::dataTableOutput("tbl_master"))
        )
      ),

      # ══════════════════════════════════════════════════════════════════
      # DOWNLOADS
      # ══════════════════════════════════════════════════════════════════
      tabItem(tabName="dl",
        fluidRow(
          box(title=LABELS$dl_title,width=12,status="primary",solidHeader=TRUE,
            fluidRow(
              column(3,tags$div(class="kpi",tags$h5(icon("table")," OG categories"),
                downloadButton("dl_cats","Download CSV",class="btn-primary"))),
              column(3,tags$div(class="kpi g",tags$h5(icon("tags")," Master annotation"),
                downloadButton("dl_ann","Download TSV",class="btn-success"))),
              column(3,tags$div(class="kpi o",tags$h5(icon("atom")," CAZymes"),
                downloadButton("dl_cazy","Download CSV",class="btn-warning"))),
              column(3,tags$div(class="kpi p",tags$h5(icon("biohazard")," Effectors"),
                downloadButton("dl_eff","Download CSV",class="btn-default")))
            ),br(),
            fluidRow(
              column(3,tags$div(class="kpi t",tags$h5(icon("shield-virus")," PHI hits"),
                downloadButton("dl_phi","Download CSV",class="btn-info"))),
              column(3,tags$div(class="kpi",tags$h5(icon("tree")," CAFE annotation"),
                downloadButton("dl_cafe","Download CSV",class="btn-default"))),
              column(3,tags$div(class="kpi",tags$h5(icon("arrow-trend-up")," Selection pressure"),
                downloadButton("dl_sel","Download CSV",class="btn-default"))),
              column(3,tags$div(class="kpi",tags$h5(icon("circle-nodes")," PPI top 10% nodes"),
                downloadButton("dl_ppi","Download CSV",class="btn-default")))
            )
          )
        )
      )
    )
  )
)

# =============================================================================
# SERVER
# =============================================================================
server <- function(input, output, session) {

  # ── Status banner ────────────────────────────────────────────────────────────
  output$status_banner <- renderUI({
    chk <- c(
      "gene_category_summary.csv" = file.exists(file.path(DATA_DIR,"gene_category_summary.csv")),
      "busco_summary.csv"         = file.exists(file.path(DATA_DIR,"busco_summary.csv")),
      "cazyme_annotations.csv"    = file.exists(file.path(VIR_DIR,"cazyme_annotations.csv")),
      "effector_annotations.csv"  = file.exists(file.path(VIR_DIR,"effector_annotations.csv")),
      "cafe_expanded_contracted_annotations.csv" = file.exists(file.path(CAFE_DIR,"cafe_expanded_contracted_annotations.csv")),
      "core_positive_selection.csv" = file.exists(file.path(SEL_DIR,"core_positive_selection.csv")),
      "core_ppi_network.svg"      = file.exists(file.path(DATA_DIR,"ppi","core_ppi_network.svg"))
    )
    cls <- if(all(chk)) "ok" else "ok warn"
    msg <- if(all(chk)) "✅  All pangenomic data tracks verified" else
      paste0("⚠️  Missing paths tracked: ",paste(names(chk)[!chk],collapse=", "))
    tags$div(class=cls, msg)
  })

  # ── HOME — World map (full width) ───────────────────────────────────────────
  output$world_map <- renderLeaflet({
    leaflet(strain_meta) |>
      addProviderTiles(providers$CartoDB.Positron) |>
      setView(lng=100, lat=20, zoom=2) |>
      addCircleMarkers(lng=~Lon, lat=~Lat, radius=13, color="#fff", weight=2,
        fillColor=~color, fillOpacity=0.9,
        popup=~paste0(
          "<b style='font-size:14px;color:#1565c0;'>",Strain,"</b><br>",
          "<b>Location:</b> ",Location,"<br>",
          "<b>Country:</b> ",Country,"<br>",
          "<b>Host:</b> ",Host,"<br>",
          "<b>Assembly:</b> ",Assembly_Level,"<br>",
          "<b>Technology:</b> ",Sequencing_Technology,"<br>",
          "<b>NCBI:</b> ",NCBI_Accession,"<br>",
          "<b>Size (Mb):</b> ",Size_Mb
        ),
        label=~Strain,
        labelOptions=labelOptions(style=list("font-weight"="600"),direction="top",offset=c(0,-15)))
  })

  # ── HOME — Strain metadata table (all 10 columns) ───────────────────────────
  output$tbl_strains <- DT::renderDataTable({
    d <- strain_meta[, c("Strain","Country","Isolation_Source","Host","Location",
                         "Assembly_Level","Pub_Year","Sequencing_Technology",
                         "NCBI_Accession","WGS_Accession","Size_Mb")]
    colnames(d) <- c("Strain","Country","Isolation Source","Host","Location",
                     "Assembly Level","Year","Sequencing Technology",
                     "NCBI Accession","WGS Accession","Size (Mb)")
    datatable(d, rownames=FALSE, filter="top",
      options=list(pageLength=10, scrollX=TRUE, autoWidth=FALSE,
                   columnDefs=list(list(width='130px', targets="_all"))))
  })

  # ── BUSCO plot — X-axis title "BUSCO %" on both panels ──────────────────────
  output$plt_busco_split <- renderPlotly({
    if(nrow(busco)==0) return(plot_ly()|>layout(title="No BUSCO source file matched"))
    b_df <- copy(busco)
    setnames(b_df, names(b_df)[1], "Strains")
    b_df <- b_df[order(Strains)]

    s1 <- plot_ly(b_df,
                  y=~Strains,
                  x=~(`Complete and single-copy BUSCOs` / `Total BUSCO` * 100),
                  type='bar', orientation='h', name='Single-Copy',
                  marker=list(color='#56B4E9')) %>%
          layout(xaxis=list(range=c(98,100), title="BUSCO %"))

    s2 <- plot_ly(b_df,
                  y=~Strains,
                  x=~(`Complete and duplicated BUSCOs` / `Total BUSCO` * 100),
                  type='bar', orientation='h', name='Duplicated',
                  marker=list(color='#3282BD')) %>%
          add_trace(x=~(`Fragmented BUSCOs`  / `Total BUSCO` * 100),
                    name='Fragmented', marker=list(color='#F0E442')) %>%
          add_trace(x=~(`Missing BUSCO` / `Total BUSCO` * 100),
                    name='Missing',    marker=list(color='#E66C4C')) %>%
          layout(barmode='stack', xaxis=list(range=c(0,2), title="BUSCO %"))

    subplot(s1, s2, nrows=1, shareY=TRUE, titleX=TRUE, margin=0.04) %>%
      layout(
        legend=list(orientation="h", y=1.1),
        margin=list(l=100, r=20, t=40, b=60)
      )
  })

  # ── PANGENOME SUMMARY ────────────────────────────────────────────────────────
  output$plt_cat_donut <- renderPlotly({
    d <- cat_counts
    plot_ly(d,labels=~Category,values=~N,type="pie",hole=0.44,
      marker=list(colors=cat_pal[as.character(d$Category)],
                  line=list(color="#fff",width=2)),
      textinfo="label+value+percent",
      hovertemplate="<b>%{label}</b>: %{value} OGs (%{percent})<extra></extra>") |>
    layout(showlegend=TRUE,
           legend=list(orientation="h",x=0,y=-0.08,font=list(size=13)),
           margin=list(t=20,b=60))
  })

  output$tbl_summary <- DT::renderDataTable({
    if(nrow(summ)==0) return(datatable(data.frame(Note="Not loaded")))
    d <- as.data.frame(summ)
    d$Total <- rowSums(d[,c("Core","Shell","Cloud","Softcore")],na.rm=TRUE)
    datatable(d,rownames=FALSE,options=list(pageLength=10,scrollX=TRUE))
  })

  output$plt_strain_stack <- renderPlotly({
    if(nrow(summ_long)==0) return(plot_ly()|>layout(title="pangenome_summary.txt not loaded"))
    totals <- summ[,.(total=Core+Shell+Cloud+Softcore),by=Genome]
    ord <- totals[order(-total),Genome]
    d <- summ_long
    d$Genome   <- factor(d$Genome,  levels=ord)
    d$Category <- factor(d$Category,levels=c("Cloud","Shell","Softcore","Core"))
    plot_ly(d,x=~Genome,y=~Count,color=~Category,
      colors=cat_pal[c("Cloud","Shell","Softcore","Core")],type="bar",
      hovertemplate="<b>%{x}</b><br>%{fullData.name}: %{y}<extra></extra>") |>
    layout(barmode="stack",xaxis=list(title="",tickangle=-40),
           yaxis=list(title="Number of genes",tickformat="d"),
           legend=list(orientation="h",x=0,y=1.08))
  })

  # ── ORTHOGROUP EXPLORER ──────────────────────────────────────────────────────
  og_tbl_data <- eventReactive(input$og_go,{
    if(nrow(cats)==0) return(data.table())
    d <- copy(cats)
    if(nchar(trimws(input$og_q))>0) d <- d[grepl(input$og_q,Orthogroup,ignore.case=TRUE)]
    if(input$og_cat!="All") d <- d[Category==input$og_cat]
    if(nrow(ann_master)>0){
      ac <- intersect(c("Orthogroup","Category","EGG_Desc","NR_Desc","IPR_Desc","SP_Hit"),names(ann_master))
      d <- merge(d,ann_master[,..ac],by=c("Orthogroup","Category"),all.x=TRUE)
    }
    d
  },ignoreNULL=FALSE)

  output$tbl_og <- DT::renderDataTable({
    datatable(og_tbl_data(),selection="single",rownames=FALSE,
      options=list(pageLength=12,scrollX=TRUE,lengthMenu=c(10,25,50,100,500)))
  })

  sel_og <- reactive({
    s <- input$tbl_og_rows_selected
    if(!length(s)) return(NULL); og_tbl_data()[s]
  })

  output$og_detail_ui <- renderUI({
    og <- sel_og(); if(is.null(og)) return(NULL)
    cc <- switch(og$Category[1],Core="pill-core",Softcore="pill-soft",Shell="pill-shell",Cloud="pill-cloud","")
    fluidRow(box(
      title=tags$span("OG: ",tags$b(og$Orthogroup[1])," ",
                      tags$span(class=paste("ann-pill",cc),og$Category[1])),
      width=12,status="success",solidHeader=TRUE,
      fluidRow(
        column(5,tags$h5(icon("dna")," Gene members"),DT::dataTableOutput("tbl_og_members")),
        column(7,tags$h5(icon("tags")," Annotation"),uiOutput("og_ann_card"),hr(),
          tags$h5(icon("dna")," Complete Sequences (Full Sequence Trace Output View)"),
          verbatimTextOutput("seq_out",placeholder=TRUE))
      )
    ))
  })

  output$seq_out <- renderText({
    og <- sel_og(); if(is.null(og)) return("Select an orthogroup to view sequences.")
    fp <- file.path(SEQ_DIR,"all_orthogroup_sequences.fa")
    if(!file.exists(fp)) return(paste0("⚠️  Sequence file path tracking missed: ",fp))
    withProgress(message=paste0("Parsing complete FASTA tracks for ",og$Orthogroup[1],"…"),{
      all_lines <- readLines(fp, warn=FALSE)
      og_pattern <- paste0("^>",og$Orthogroup[1],"\\|")
      seq_indices <- grep(og_pattern, all_lines)
      if(length(seq_indices)==0) return(paste0("No matching FASTA records located for ",og$Orthogroup[1]))
      output_lines <- character()
      for(i in seq_indices) {
        output_lines <- c(output_lines, all_lines[i])
        if(i < length(all_lines)) {
          j <- i+1
          while(j<=length(all_lines) && !grepl("^>",all_lines[j])){ output_lines <- c(output_lines,all_lines[j]); j <- j+1 }
        }
      }
    })
    paste(output_lines,collapse="\n")
  })

  output$tbl_og_members <- DT::renderDataTable({
    og <- sel_og(); if(is.null(og)) return(NULL)
    f <- file.path(DATA_DIR,"orthogroup_summary.tsv")
    if(!file.exists(f)) return(datatable(data.frame(Note="TSV map tracking error")))
    withProgress(message="Extracting strain locus tags…",{
      dt_full <- fread(f,sep="\t",header=TRUE,fill=TRUE,data.table=TRUE)
      row <- dt_full[Orthogroup==og$Orthogroup[1]]
    })
    if(nrow(row)==0) return(datatable(data.frame(Note="OG not found")))
    sc <- setdiff(names(row),"Orthogroup")
    long <- data.frame(Strain=sub("_renamed$","",sc),
      Genes=sapply(sc,function(c) trimws(as.character(row[[c]]))),stringsAsFactors=FALSE)
    long$N_Genes <- sapply(long$Genes,function(g) if(nchar(g)==0) 0L else length(strsplit(g,",")[[1]]))
    datatable(long[long$N_Genes>0,c("Strain","N_Genes","Genes")],rownames=FALSE,
      options=list(pageLength=8,scrollX=TRUE))
  })

  output$og_ann_card <- renderUI({
    og <- sel_og(); if(is.null(og)||nrow(ann_master)==0) return(NULL)
    r <- ann_master[Orthogroup==og$Orthogroup[1]]
    if(nrow(r)==0) return(tags$p(style="color:#888;","No annotation."))
    gf <- function(col) if(col %in% names(r) && !is.na(r[[col]][1])) r[[col]][1] else NULL
    items <- list()
    if(!is.null(gf("EGG_Desc"))) items[[length(items)+1]] <- tags$p(tags$b("EggNOG: "),gf("EGG_Desc"))
    if(!is.null(gf("NR_Desc")))  items[[length(items)+1]] <- tags$p(tags$b("NRDB: "),gf("NR_Desc"))
    if(!is.null(gf("IPR_Desc"))) items[[length(items)+1]] <- tags$p(tags$b("InterPro: "),gf("IPR_Desc"))
    if(!is.null(gf("SP_Hit")))   items[[length(items)+1]] <- tags$p(tags$b("SwissProt: "),gf("SP_Hit"))
    if(!is.null(gf("GOs")))      items[[length(items)+1]] <- tags$p(tags$b("GO: "),gf("GOs"))
    if(!length(items)) return(tags$p(style="color:#888;","No annotation."))
    do.call(tagList,items)
  })

  # ── PRESENCE / ABSENCE ───────────────────────────────────────────────────────
  PAGE_SIZE_PA <- 100L
  pa_page <- reactiveVal(1L)
  pa_mode <- reactiveVal("page")
  observeEvent(input$pa_cat,{ pa_page(1L); pa_mode("page") })
  observeEvent(input$pa_prev,{ pa_page(max(1L,pa_page()-1L)) })
  observeEvent(input$pa_next,{
    d <- pa_filtered_all()
    if(!is.null(d)) pa_page(min(ceiling(nrow(d)/PAGE_SIZE_PA),pa_page()+1L))
  })
  observeEvent(input$pa_search_btn,{ pa_mode("search"); pa_page(1L) })

  pa_all_data <- reactive({
    if(nrow(pa_pre)==0) return(NULL)
    d <- copy(pa_pre); setnames(d,names(d)[1],"Gene_ID"); d
  })
  pa_filtered_all <- reactive({
    d <- pa_all_data(); if(is.null(d)) return(NULL)
    if(input$pa_cat!="All" && "Category" %in% names(d)) d <- d[Category==input$pa_cat]
    d
  })

  output$pa_page_info <- renderUI({
    d <- pa_filtered_all(); if(is.null(d)) return(NULL)
    n_pages <- ceiling(nrow(d)/PAGE_SIZE_PA)
    tags$div(class="heat-info",
      sprintf("Page %d / %d · %s total blocks mapped",pa_page(),n_pages,fmt(nrow(d))))
  })
  output$pa_box_title <- renderUI({
    if(pa_mode()=="search")
      tags$span("Presence / absence Matrix trace — query: ",tags$b(input$pa_search))
    else {
      d <- pa_filtered_all()
      np <- if(!is.null(d)) ceiling(nrow(d)/PAGE_SIZE_PA) else 0
      tags$span(sprintf("Presence/Absence of Orthogroup — Category: %s (Block %d/%d)",input$pa_cat,pa_page(),np))
    }
  })
  output$pa_search_result <- renderUI({
    req(pa_mode()=="search",nchar(trimws(input$pa_search))>0)
    d <- pa_all_data(); if(is.null(d)) return(NULL)
    r <- d[grepl(input$pa_search,Gene_ID,ignore.case=TRUE)]
    if(nrow(r)==0) return(tags$div(class="warn",style="padding:8px;border-radius:5px;","Locus block absent."))
    r1 <- r[1]
    sc <- setdiff(names(r1),c("Gene_ID","Category","presence_count"))
    pres <- sc[sapply(sc,function(c) !is.na(r1[[c]])&&suppressWarnings(as.integer(r1[[c]]))==1L)]
    abs_ <- setdiff(sc,pres)
    tags$div(class="heat-info",style="margin-bottom:8px;",
      tags$b(r1$Gene_ID[1])," · ",tags$b(r1$Category[1]),br(),
      tags$b("Present across strains: "),paste(pres,collapse=", "),br(),
      tags$b("Absent strains: "),paste(abs_,collapse=", "))
  })

  output$plt_pa <- renderPlotly({
    if(pa_mode()=="search"&&nchar(trimws(input$pa_search))>0){
      d <- pa_all_data(); if(is.null(d)) return(plot_ly())
      d <- d[grepl(input$pa_search,Gene_ID,ignore.case=TRUE)]
      if(nrow(d)==0) return(plot_ly()|>layout(title="Locus record not matched"))
    } else {
      d <- pa_filtered_all(); if(is.null(d)) return(plot_ly())
      pg <- pa_page()
      d <- d[(((pg-1L)*PAGE_SIZE_PA)+1L):min(pg*PAGE_SIZE_PA,nrow(d))]
    }
    if(is.null(d)||nrow(d)==0) return(plot_ly())
    s_cols <- setdiff(names(d),c("Gene_ID","Category","presence_count"))
    mat_list <- lapply(s_cols,function(c) suppressWarnings(as.numeric(d[[c]])))
    mat <- do.call(cbind,mat_list); colnames(mat) <- s_cols
    plot_ly(x=s_cols,y=as.character(d$Gene_ID),z=mat,type="heatmap",
      colorscale=list(c(0,"#ffffff"),c(1,"#1565c0")),showscale=TRUE,
      colorbar=list(title="",tickvals=c(0,1),ticktext=c("Absent","Present"),len=0.25),
      hovertemplate="Strain: <b>%{x}</b><br>OG: <b>%{y}</b><extra></extra>") |>
    layout(xaxis=list(title="",tickangle=-45,tickfont=list(size=10)),
           yaxis=list(title="",tickfont=list(size=9),dtick=1,autorange="reversed"),
           margin=list(l=160,b=130,t=10))
  })

  # ── ANNOTATION HUB ───────────────────────────────────────────────────────────
  output$plt_ann_coverage <- renderPlotly({
    if(nrow(ann_summ)>0){
      d <- copy(ann_summ); setnames(d,names(d)[1],"Category")
      d$Category <- stringr::str_to_title(trimws(as.character(d$Category)))
      d$Category <- factor(d$Category,levels=c("Core","Softcore","Shell","Cloud"))
      d <- d[order(d$Category),]
      for(i in seq_along(names(d))){ lw <- tolower(trimws(names(d)[i]))
        if(grepl("total",lw))     setnames(d,names(d)[i],"Total orthogroups")
        if(grepl("nrdb",lw))      setnames(d,names(d)[i],"NRDB")
        if(grepl("eggnog",lw))    setnames(d,names(d)[i],"EggNOG")
        if(grepl("interpro",lw))  setnames(d,names(d)[i],"InterProScan")
        if(grepl("swissprot",lw)) setnames(d,names(d)[i],"SwissProt")
      }
      cn2 <- setdiff(names(d),"Category")
      dl <- melt(d,id.vars="Category",measure.vars=cn2,variable.name="Source",value.name="Count")
      dl$Count <- suppressWarnings(as.numeric(dl$Count))
      src_col <- c("Total orthogroups"="#1565C0","NRDB"="#2E7D32","EggNOG"="#E65100",
                   "InterProScan"="#6A1B9A","SwissProt"="#00695C")
      plot_ly(dl,x=~Category,y=~Count,color=~Source,
        colors=src_col[intersect(names(src_col),unique(as.character(dl$Source)))],type="bar",
        hovertemplate="<b>%{x}</b> — %{fullData.name}: %{y}<extra></extra>") |>
      layout(barmode="group",xaxis=list(title=""),
             yaxis=list(title="Orthogroups",tickformat="d"),
             legend=list(orientation="h",x=0,y=1.12))
    } else plot_ly()|>layout(title="Source matrix unread")
  })

  ann_filt <- eventReactive(input$ann_go,{
    d <- copy(ann_master)
    if(input$ann_cat_f!="All") d <- d[Category==input$ann_cat_f]
    if(nchar(trimws(input$ann_q))>0){
      q <- trimws(input$ann_q)
      cs <- intersect(c("Orthogroup","EGG_Desc","EGG_Name","NR_Desc","IPR_Desc","SP_Hit"),names(d))
      hits <- apply(d[,cs,with=FALSE],1,function(r) any(grepl(q,r,ignore.case=TRUE)))
      d <- d[hits]
    }
    d
  },ignoreNULL=FALSE)

  output$tbl_ann <- DT::renderDataTable({
    d <- ann_filt()
    sc <- intersect(c("Orthogroup","Category","COG","EGG_Desc","EGG_Name","NR_Desc",
                      "NR_Species","IPR_Acc","IPR_Desc","SP_Hit","SP_pident","GOs","KEGG","PFAMs"),names(d))
    datatable(d[,sc,with=FALSE],rownames=FALSE,filter="top",
      options=list(pageLength=12,scrollX=TRUE))
  })

  # ── CAZymes ───────────────────────────────────────────────────────────────────
  output$plt_cazy_strain <- renderPlotly({
    if(nrow(cazy_og)==0) return(plot_ly()|>layout(title="Not loaded"))
    d <- melt(cazy_og,id.vars=c("Strain","Total_CAZymes"),
              measure.vars=intersect(c("Core","Softcore","Shell","Cloud"),names(cazy_og)),
              variable.name="Category",value.name="Count")
    d$Strain   <- factor(d$Strain,  levels=cazy_og[order(-Total_CAZymes),Strain])
    d$Category <- factor(d$Category,levels=c("Cloud","Shell","Softcore","Core"))
    plot_ly(d,x=~Strain,y=~Count,color=~Category,
      colors=cat_pal[c("Cloud","Shell","Softcore","Core")],type="bar",
      hovertemplate="<b>%{x}</b><br>%{fullData.name}: %{y}<extra></extra>") |>
    layout(barmode="stack",xaxis=list(title="",tickangle=-40),
           yaxis=list(title="CAZyme genes",tickformat="d"),
           legend=list(orientation="h",x=0,y=1.08))
  })

  CAZY_PAGE_CUSTOM <- 30L
  cazy_heat_page <- reactiveVal(1L)
  cazy_heat_filtered <- reactiveVal(NULL)

  observeEvent(input$cazy_heat_go,{
    cazy_heat_page(1L)
    if(nrow(cazy_sub)==0){ cazy_heat_filtered(NULL); return() }
    d <- copy(cazy_sub)
    q <- trimws(input$cazy_heat_q)
    if(nchar(q)>0) d <- d[grepl(q,Subclass,ignore.case=TRUE)]
    cazy_heat_filtered(d[order(-Total)])
  },ignoreNULL=FALSE)

  cazy_heat_data_rv <- reactive({
    d <- cazy_heat_filtered()
    if(is.null(d)){ if(nrow(cazy_sub)>0) return(cazy_sub[order(-Total)]) else return(NULL) }
    d
  })

  observeEvent(input$cazy_heat_prev,{ cazy_heat_page(max(1L,cazy_heat_page()-1L)) })
  observeEvent(input$cazy_heat_next,{
    d <- cazy_heat_data_rv()
    if(!is.null(d)) cazy_heat_page(min(ceiling(nrow(d)/CAZY_PAGE_CUSTOM),cazy_heat_page()+1L))
  })

  output$cazy_heat_page_info <- renderUI({
    d <- cazy_heat_data_rv(); if(is.null(d)) return(NULL)
    np <- ceiling(nrow(d)/CAZY_PAGE_CUSTOM)
    tags$div(class="heat-info",sprintf("Page %d / %d · %d total families mapped",cazy_heat_page(),np,nrow(d)))
  })

  output$plt_cazy_heat <- renderPlotly({
    d <- cazy_heat_data_rv()
    if(is.null(d)||nrow(d)==0) return(plot_ly()|>layout(title="No structures parsed"))
    pg <- cazy_heat_page()
    i1 <- (pg-1L)*CAZY_PAGE_CUSTOM+1L; i2 <- min(pg*CAZY_PAGE_CUSTOM,nrow(d))
    d_page <- d[i1:i2]
    s_cols <- setdiff(names(d_page),c("Subclass","Total"))
    mat <- as.matrix(d_page[,s_cols,with=FALSE]); mode(mat) <- "numeric"
    plot_ly(x=s_cols,y=as.character(d_page$Subclass),z=mat,type="heatmap",
      colorscale=list(c(0,"#FFFFFF"),c(0.5,"#FFD700"),c(1,"#B22222")),
      showscale=TRUE,colorbar=list(title="Count"),zmin=0,
      hovertemplate="Strain: <b>%{x}</b><br>Family: <b>%{y}</b><br>Count: %{z}<extra></extra>") |>
    layout(xaxis=list(title="",tickangle=-45,tickfont=list(size=11)),
           yaxis=list(title="",tickfont=list(size=10),dtick=1,autorange="reversed"),
           margin=list(l=140,b=120,t=50))
  })

  output$plt_cazy_bar <- renderPlotly({
    if(nrow(cazy_sub)==0) return(plot_ly()|>layout(title="Not loaded"))
    d <- cazy_sub[order(-Total)]
    n <- nrow(d)
    plasma_cols <- c("#0D0887FF","#1E1287FF","#3D0A7CFF","#5B0C75FF","#78096BFF",
                     "#950A6EFF","#B2186FFF","#CB236EFF","#E23E6AFF","#F85C62FF",
                     "#FEB13CFF","#FEDA40FF")
    cls <- colorRampPalette(plasma_cols)(n)
    plot_ly(d,x=~reorder(Subclass,-Total),y=~Total,type="bar",
      marker=list(color=cls,line=list(color="#ccc",width=0.3)),
      hovertemplate="<b>%{x}</b>: %{y} total<extra></extra>") |>
    layout(xaxis=list(title="CAZyme subclass",tickangle=-70,tickfont=list(size=11)),
           yaxis=list(title="Total count (all strains)",tickformat="d"),margin=list(b=180))
  })

  cazy_filt <- reactive({
    d <- copy(cazy)
    if("Strain" %in% names(d) && input$cazy_f_strain!="All") d <- d[Strain==input$cazy_f_strain]
    if("Pangenome_Category" %in% names(d) && input$cazy_f_cat!="All") d <- d[Pangenome_Category==input$cazy_f_cat]
    if("CAZyme_Subclass" %in% names(d) && input$cazy_f_sub!="All") d <- d[CAZyme_Subclass==input$cazy_f_sub]
    d
  })
  output$tbl_cazy <- DT::renderDataTable({
    datatable(cazy_filt(),rownames=FALSE,filter="top",options=list(pageLength=12,scrollX=TRUE))
  })

  # ── EFFECTORS ─────────────────────────────────────────────────────────────────
  output$plt_eff_type <- renderPlotly({
    if(nrow(eff_summ)==0) return(plot_ly()|>layout(title="Not loaded"))
    tc <- intersect(c("Apoplastic effector","Cytoplasmic effector","Dual-effector"),names(eff_summ))
    d <- melt(eff_summ,id.vars=c("Strain","total_effector"),measure.vars=tc,
              variable.name="Type",value.name="Count")
    d$Strain <- factor(d$Strain,levels=eff_summ[order(-total_effector),Strain])
    plot_ly(d,x=~Strain,y=~Count,color=~Type,
      colors=c("#1565c0","#e65100","#6a1b9a"),type="bar",
      hovertemplate="<b>%{x}</b><br>%{fullData.name}: %{y}<extra></extra>") |>
    layout(barmode="stack",xaxis=list(title="",tickangle=-40),
           yaxis=list(title="Effectors",tickformat="d"),legend=list(orientation="h",x=0,y=1.08))
  })

  output$plt_eff_cat <- renderPlotly({
    if(nrow(eff)==0) return(plot_ly()|>layout(title="Not loaded"))
    d <- copy(eff)
    if("Pangenome_Category" %in% names(d)){
      d[,Pangenome_Category:=stringr::str_to_title(trimws(Pangenome_Category))]
      d <- d[,.N,by=Pangenome_Category]
      plot_ly(d,labels=~Pangenome_Category,values=~N,type="pie",hole=0.45,
        marker=list(colors=cat_pal[as.character(d$Pangenome_Category)],
                    line=list(color="#fff",width=2)),
        hovertemplate="<b>%{label}</b>: %{value} (%{percent})<extra></extra>") |>
      layout(showlegend=TRUE,margin=list(t=10,b=10))
    } else plot_ly()
  })

  eff_filt <- reactive({
    d <- copy(eff)
    if("Contributing_Strains" %in% names(d) && input$eff_f_strain!="All")
      d <- d[grepl(input$eff_f_strain,Contributing_Strains,ignore.case=TRUE)]
    if("Pangenome_Category" %in% names(d) && input$eff_f_cat!="All")
      d <- d[stringr::str_to_title(trimws(Pangenome_Category))==input$eff_f_cat]
    if("Effector_Type" %in% names(d) && input$eff_f_type!="All")
      d <- d[Effector_Type==input$eff_f_type]
    d
  })
  output$tbl_eff <- DT::renderDataTable({
    datatable(eff_filt(),rownames=FALSE,filter="top",options=list(pageLength=12,scrollX=TRUE))
  })

  # ── PHI DATABASE ──────────────────────────────────────────────────────────────
  output$plt_phi_strain <- renderPlotly({
    if(nrow(phi_strain)==0) return(plot_ly()|>layout(title="Not loaded"))
    d <- copy(phi_strain); setnames(d,names(d)[1],"Strain")
    phen_cols <- setdiff(names(d),c("Strain","Total_PHI_Hits"))
    if("Total_PHI_Hits" %in% names(d)) d <- d[order(-Total_PHI_Hits)]
    d$Strain <- factor(d$Strain,levels=d$Strain)
    dl <- melt(d,id.vars="Strain",measure.vars=phen_cols,variable.name="Phenotype",value.name="Count")
    dl$Count <- suppressWarnings(as.numeric(dl$Count))
    n <- length(phen_cols)
    cols <- colorRampPalette(brewer.pal(min(n,11),"Spectral"))(n)
    plot_ly(dl,x=~Strain,y=~Count,color=~Phenotype,colors=setNames(cols,phen_cols),type="bar",
      hovertemplate="<b>%{x}</b><br>%{fullData.name}: %{y}<extra></extra>") |>
    layout(barmode="stack",xaxis=list(title="",tickangle=-40),
           yaxis=list(title="PHI hits",tickformat="d"),
           legend=list(orientation="v",x=1.01,y=1,font=list(size=12)),
           margin=list(r=360,b=80))
  })

  output$tbl_phi <- DT::renderDataTable({
    if(nrow(phi_db)==0) return(datatable(data.frame(Note="Not loaded")))
    sc <- if(all(c("OG","Gene_ID","PHI_Gene","Organism","Phenotype") %in% names(phi_db)))
      c("OG","Gene_ID","PHI_Gene","Organism","Phenotype",
        intersect(c("pident","length","evalue","bitscore"),names(phi_db)))
    else names(phi_db)
    datatable(phi_db[,sc,with=FALSE],rownames=FALSE,filter="top",
      options=list(pageLength=12,scrollX=TRUE))
  })

  # ── CAFE ──────────────────────────────────────────────────────────────────────
  output$plt_cafe_ec <- renderPlotly({
    if(nrow(cafe_sum)==0) return(plot_ly()|>layout(title="Not loaded"))
    d <- melt(cafe_sum,id.vars="Strain",variable.name="Type",value.name="Count")
    d$Strain <- factor(d$Strain,levels=cafe_sum[order(-(Expansions+Contractions)),Strain])
    plot_ly(d,x=~Strain,y=~Count,color=~Type,
      colors=c(Expansions="#2E7D32",Contractions="#C62828"),type="bar",
      hovertemplate="<b>%{x}</b><br>%{fullData.name}: %{y}<extra></extra>") |>
    layout(barmode="group",xaxis=list(title="",tickangle=-40),
           yaxis=list(title="Gene families",tickformat="d"),legend=list(orientation="h",x=0,y=1.08))
  })

  output$plt_cafe_cat <- renderPlotly({
    if(nrow(cafe_ann)==0) return(plot_ly()|>layout(title="Not loaded"))
    d <- cafe_ann[,.N,by=Category]
    plot_ly(d,labels=~Category,values=~N,type="pie",hole=0.44,
      marker=list(colors=c("#2E7D32","#C62828","#FF9800","#1565C0")[seq_len(nrow(d))],
                  line=list(color="#fff",width=2)),
      hovertemplate="<b>%{label}</b>: %{value} OGs (%{percent})<extra></extra>") |>
    layout(showlegend=TRUE,margin=list(t=10,b=10))
  })

  CAFE_PAGE_CUSTOM <- 25L
  cafe_heat_page <- reactiveVal(1L)
  observeEvent(input$cafe_heat_prev,{ cafe_heat_page(max(1L,cafe_heat_page()-1L)) })
  observeEvent(input$cafe_heat_next,{
    if(nrow(cafe_heat)>0) cafe_heat_page(min(ceiling(nrow(cafe_heat)/CAFE_PAGE_CUSTOM),cafe_heat_page()+1L))
  })

  output$cafe_heat_page_info <- renderUI({
    if(nrow(cafe_heat)==0) return(NULL)
    np <- ceiling(nrow(cafe_heat)/CAFE_PAGE_CUSTOM)
    tags$div(class="heat-info",sprintf("Page %d / %d · %d parsed matrix sets",cafe_heat_page(),np,nrow(cafe_heat)))
  })

  output$plt_cafe_heat <- renderPlotly({
    if(nrow(cafe_heat)==0) return(plot_ly()|>layout(title="No tracks read"))
    d <- copy(cafe_heat); setnames(d,names(d)[1],"FamilyID")
    pg <- cafe_heat_page()
    i1 <- (pg-1L)*CAFE_PAGE_CUSTOM+1L; i2 <- min(pg*CAFE_PAGE_CUSTOM,nrow(d))
    d_page <- d[i1:i2]
    s_cols <- setdiff(names(d_page),"FamilyID")
    mat <- as.matrix(d_page[,s_cols,with=FALSE]); mode(mat) <- "numeric"; rownames(mat) <- NULL
    plot_ly(x=s_cols,y=as.character(d_page$FamilyID),z=mat,type="heatmap",
      colorscale=list(c(0,"#FFFFFF"),c(0.5,"#FFD700"),c(1,"#B22222")),
      showscale=TRUE,colorbar=list(title="Gene copies"),zmin=0,
      hovertemplate="Strain: <b>%{x}</b><br>OG: <b>%{y}</b><br>Copies: %{z}<extra></extra>") |>
    layout(xaxis=list(title="",tickangle=-45,tickfont=list(size=11)),
           yaxis=list(title="",tickfont=list(size=10),dtick=1,autorange="reversed"),
           margin=list(l=140,b=120,t=50))
  })

  output$tbl_cafe <- DT::renderDataTable({
    if(nrow(cafe_ann)==0) return(datatable(data.frame(Note="Not loaded")))
    datatable(cafe_ann,rownames=FALSE,filter="top",options=list(pageLength=12,scrollX=TRUE))
  })

  # ── SELECTION PRESSURE ────────────────────────────────────────────────────────
  output$plt_sel_summary <- renderPlotly({
    if(nrow(sel_core)==0) return(plot_ly()|>layout(title="Not loaded"))
    d <- sel_core[,.N,by=Selection]
    plot_ly(d,labels=~Selection,values=~N,type="pie",hole=0.44,
      marker=list(colors=c("#C62828","#2E7D32","#1565C0","#E65100")[seq_len(nrow(d))],
                  line=list(color="#fff",width=2)),
      hovertemplate="<b>%{label}</b>: %{value} OGs (%{percent})<extra></extra>") |>
    layout(showlegend=TRUE,title=list(text="Selection profile tracks",x=0.5))
  })

  output$plt_sel_methods <- renderPlotly({
    method_data <- data.frame(
      Method=c("FEL(+)","FEL(-)","MEME","BUSTED","aBSREL"),
      Count =c(619,6554,4542,4408,6927))
    method_colors <- c("#2E7D32","#E65100","#1565C0","#6A1B9A","#FF9800")
    plot_ly(method_data,x=~reorder(Method,-Count),y=~Count,type="bar",
      marker=list(color=method_colors),
      hovertemplate="<b>%{x}</b>: %{y} orthogroups<extra></extra>") |>
    layout(xaxis=list(title="Selection method"),
           yaxis=list(title="Number of orthogroups",tickformat="d"),margin=list(b=60))
  })

  output$tbl_sel_strong <- DT::renderDataTable({
    if(nrow(sel_strong)==0) return(datatable(data.frame(Note="Not loaded")))
    datatable(sel_strong,rownames=FALSE,filter="top",options=list(pageLength=10,scrollX=TRUE))
  })
  output$tbl_sel_all <- DT::renderDataTable({
    if(nrow(sel_core)==0) return(datatable(data.frame(Note="Not loaded")))
    datatable(sel_core,rownames=FALSE,filter="top",options=list(pageLength=10,scrollX=TRUE))
  })

  # ── CORE PPI NETWORK ─────────────────────────────────────────────────────────
  ppi_state <- reactiveVal("top10")
  ppi_search_result <- reactiveVal(NULL)

  observeEvent(input$ppi_reset_btn,{ ppi_state("top10"); ppi_search_result(NULL) })

  observeEvent(input$ppi_search_btn,{
    og_q <- trimws(input$ppi_search_og)
    if(nchar(og_q)==0){ ppi_state("top10"); return() }
    if(nrow(ppi_ann_full)==0){ showNotification("Master annotations unread",type="warning"); return() }
    an_col <- names(ppi_ann_full)[1]
    node_row <- ppi_ann_full[grepl(og_q,get(an_col),ignore.case=TRUE)]
    if(nrow(node_row)==0){ showNotification(paste("ID not tracked:",og_q),type="warning"); return() }
    target_id <- as.character(node_row[[an_col]][1])
    full_edge_f <- file.path(PPI_CORE,"core_ppi_interactions_score400.tsv")
    if(!file.exists(full_edge_f)){ showNotification("Sub-interaction matrix unlinked",type="warning"); return() }
    withProgress(message=paste0("Parsing connections for ",target_id,"…"),{
      all_e <- fread(full_edge_f,sep="\t",header=TRUE,data.table=FALSE)
      names(all_e)[1:2] <- c("from","to")
      sub_e <- all_e[all_e$from==target_id | all_e$to==target_id, ]
    })
    if(nrow(sub_e)==0){ showNotification("No interaction node vectors unmasked",type="warning"); return() }
    ppi_search_result(list(target=target_id,edges=sub_e,ann=ppi_ann_full))
    ppi_state("search")
  })

  output$ppi_net_title <- renderUI({
    if(ppi_state()=="search" && !is.null(ppi_search_result()))
      tags$span("Interactive Neighborhood Subnetwork Map for: ",tags$b(ppi_search_result()$target))
    else tags$span(LABELS$ppi_network_title)
  })

  output$ui_ppi_canvas <- renderUI({
    if(ppi_state()=="search" && !is.null(ppi_search_result())){
      visNetworkOutput("plt_ppi_interactive",height="620px")
    } else {
      tags$div(style="text-align:center;background:#ffffff;padding:10px;border-radius:8px;",
        tags$img(src="ppi_mesh/core_ppi_network.svg",
                 style="max-width:100%;max-height:600px;object-fit:contain;"))
    }
  })

  output$plt_ppi_interactive <- visNetwork::renderVisNetwork({
    req(ppi_state()=="search",!is.null(ppi_search_result()))
    res   <- ppi_search_result()
    sub_e <- res$edges
    names(sub_e)[1:2] <- c("from","to")
    if(ncol(sub_e)>=3) sub_e$value <- suppressWarnings(as.numeric(sub_e[[3]]))
    node_ids <- unique(c(sub_e$from,sub_e$to))
    an_col   <- names(res$ann)[1]
    nd <- data.frame(id=node_ids,label=node_ids,stringsAsFactors=FALSE)
    ann_sub <- as.data.frame(res$ann)[as.character(res$ann[[an_col]]) %in% node_ids,]
    if(nrow(ann_sub)>0){ names(ann_sub)[1] <- "id"; nd <- merge(nd,ann_sub,by="id",all.x=TRUE) }
    ann_v <- if("Annotation" %in% names(nd)) as.character(nd$Annotation) else rep("Other proteins",nrow(nd))
    ann_v[is.na(ann_v)] <- "Other proteins"
    node_colors <- rep(ann_type_pal["Other proteins"],nrow(nd))
    for(idx in seq_len(nrow(nd))){
      ts <- trimws(ann_v[idx])
      if(ts %in% names(ann_type_pal)) node_colors[idx] <- ann_type_pal[ts]
    }
    nd$color <- node_colors
    nd$color[nd$id==res$target] <- "#FFD600"
    nd$value <- 22; nd$value[nd$id==res$target] <- 55
    nd$title <- paste0("<b>",nd$id,"</b><br>Track Profile: ",ann_v)
    visNetwork::visNetwork(nd,sub_e) |>
      visNetwork::visOptions(highlightNearest=list(enabled=TRUE,degree=1),nodesIdSelection=FALSE) |>
      visNetwork::visPhysics(stabilization=list(iterations=180),enabled=TRUE) |>
      visNetwork::visEdges(smooth=FALSE,color=list(color="#9e9e9e",opacity=0.6),width=1.2) |>
      visNetwork::visNodes(font=list(size=12,face="sans-serif")) |>
      visNetwork::visInteraction(navigationButtons=TRUE,tooltipDelay=100)
  })

  output$ppi_legend_ui <- renderUI({
    items <- lapply(names(ann_type_pal),function(nm){
      tags$div(style="margin-bottom:3px;font-size:11px;",
        tags$span(class="legend-dot",style=paste0("background:",ann_type_pal[nm],";")),nm)
    })
    do.call(tagList,items)
  })

  # PPI table — ALL columns from orthogroup_annotations.csv
  ppi_tbl_rv <- eventReactive(input$ppi_tbl_go,{
    d <- if(nrow(ppi_ann_full)>0) copy(ppi_ann_full) else data.table()
    if(nrow(d)==0) return(d)
    q <- trimws(input$ppi_tbl_q)
    if(nchar(q)>0){
      char_cols <- names(d)[sapply(d,function(col) is.character(col)||is.factor(col))]
      if(length(char_cols)>0){
        hits <- apply(as.data.frame(d[,char_cols,with=FALSE]),1,
                      function(r) any(grepl(q,r,ignore.case=TRUE)))
        d <- d[hits]
      }
    }
    d
  },ignoreNULL=FALSE)

  output$tbl_ppi_ann <- DT::renderDataTable({
    d <- ppi_tbl_rv()
    if(nrow(d)==0) d <- if(nrow(ppi_ann_full)>0) copy(ppi_ann_full) else
      data.table(Status="No dataset matched")
    datatable(d, rownames=FALSE, filter="top",
      options=list(pageLength=12, scrollX=TRUE, scrollY="500px", autoWidth=FALSE,
                   columnDefs=list(list(width='120px',targets="_all"))))
  })

  # ── PANGENOME OVERVIEW ────────────────────────────────────────────────────────
  output$plt_func_cat_dist <- renderPlotly({
    func_matrix <- data.table(
      Category=c("Core","Softcore","Shell","Cloud"),
      CAFE    =c(435,0,0,0),
      CAZymes =c(290,20,180,108),
      Effector=c(95,25,210,153),
      PHI     =c(1710,130,850,353),
      Selection=c(501,0,0,0),
      PPI     =c(645,0,0,0))
    dl <- melt(func_matrix,id.vars="Category",
               measure.vars=c("CAFE","CAZymes","Effector","PHI","Selection","PPI"),
               variable.name="Functional_Class",value.name="Count")
    dl$Category <- factor(dl$Category,levels=c("Core","Softcore","Shell","Cloud"))
    levels(dl$Functional_Class)[levels(dl$Functional_Class)=="PPI"] <- "Top 10% PPI Nodes"
    tick_values <- c(0,50,100,200,300,400,500,1000,1500,2000)
    plot_ly(dl,y=~Category,x=~Count,color=~Functional_Class,
      colors=c(CAFE="#FF6F00",CAZymes="#1B5E20",Effector="#1565C0",PHI="#6A1B9A",
               Selection="#E91E63",`Top 10% PPI Nodes`="#00695C"),
      type="bar",orientation="h",
      hovertemplate="Category: <b>%{y}</b><br>%{fullData.name}: %{customdata}<extra></extra>",
      customdata=~Count) |>
    layout(barmode="group",
           xaxis=list(title="Orthogroup Density Matrix Scale Range",type="linear",
                      tickvals=tick_values,ticktext=as.character(tick_values)),
           yaxis=list(title=""),
           legend=list(orientation="h",x=0,y=1.12,font=list(size=10)))
  })

  master_data <- eventReactive(input$master_go,{
    if(nrow(cats)==0) return(data.table())
    d <- copy(cats)
    d[,CAFE            := Orthogroup %in% cafe_ogs]
    d[,CAZyme          := Orthogroup %in% cazy_ogs]
    d[,Effector        := Orthogroup %in% eff_ogs]
    d[,PHI             := Orthogroup %in% phi_ogs]
    d[,PositiveSelection:=Orthogroup %in% sel_ogs]
    d[,PPI_Top10       := Orthogroup %in% top10_ogs]
    if(nrow(ann_master)>0){
      ac <- intersect(c("Orthogroup","EGG_Desc","NR_Desc","IPR_Desc","SP_Hit"),names(ann_master))
      if(length(ac)>1) d <- merge(d,ann_master[,..ac],by="Orthogroup",all.x=TRUE)
    }
    if(input$master_cat!="All") d <- d[Category==input$master_cat]
    if(nchar(trimws(input$master_q))>0){
      q <- trimws(input$master_q)
      cs <- intersect(c("Orthogroup","EGG_Desc","NR_Desc","IPR_Desc","SP_Hit"),names(d))
      hits <- apply(d[,cs,with=FALSE],1,function(r) any(grepl(q,r,ignore.case=TRUE)))
      d <- d[hits]
    }
    if(input$master_func!="All"){
      if(input$master_func=="Multi"){
        d <- d[rowSums(d[,c("CAFE","CAZyme","Effector","PHI","PositiveSelection","PPI_Top10"),with=FALSE])>1]
      } else if(input$master_func=="None"){
        d <- d[!CAFE & !CAZyme & !Effector & !PHI & !PositiveSelection & !PPI_Top10]
      } else {
        col_map <- c(CAFE="CAFE",CAZymes="CAZyme",Effector="Effector",PHI="PHI",
                     Selection="PositiveSelection",PPI="PPI_Top10")
        col_name <- col_map[input$master_func]
        if(!is.na(col_name)) d <- d[get(col_name)==TRUE]
      }
    }
    d
  },ignoreNULL=FALSE)

  output$tbl_master <- DT::renderDataTable({
    d <- master_data()
    if(nrow(d)==0) return(datatable(data.frame(Note="No orthogroups match the filters")))
    sc <- intersect(c("Orthogroup","Category","CAFE","CAZyme","Effector","PHI",
                      "PositiveSelection","PPI_Top10","EGG_Desc","NR_Desc","IPR_Desc","SP_Hit"),names(d))
    for(col in c("CAFE","CAZyme","Effector","PHI","PositiveSelection","PPI_Top10"))
      if(col %in% names(d)) d[[col]] <- ifelse(d[[col]],"TRUE","FALSE")
    datatable(d[,sc,with=FALSE],rownames=FALSE,filter="top",options=list(pageLength=25,scrollX=TRUE))
  })

  # ── DOWNLOADS ─────────────────────────────────────────────────────────────────
  output$dl_cats <- downloadHandler("og_categories.csv",     function(f) if(nrow(cats)>0)       fwrite(cats,f))
  output$dl_ann  <- downloadHandler("annotation_master.tsv", function(f) if(nrow(ann_master)>0) fwrite(ann_master,f,sep="\t"))
  output$dl_cazy <- downloadHandler("cazymes.csv",           function(f) if(nrow(cazy)>0)       fwrite(cazy,f))
  output$dl_eff  <- downloadHandler("effectors.csv",         function(f) if(nrow(eff)>0)        fwrite(eff,f))
  output$dl_phi  <- downloadHandler("phi_hits.csv",          function(f) if(nrow(phi_db)>0)     fwrite(phi_db,f))
  output$dl_cafe <- downloadHandler("cafe_annotation.csv",   function(f) if(nrow(cafe_ann)>0)   fwrite(cafe_ann,f))
  output$dl_sel  <- downloadHandler("selection_pressure.csv",function(f) if(nrow(sel_core)>0)   fwrite(sel_core,f))
  output$dl_ppi  <- downloadHandler("ppi_top10_nodes.csv",   function(f) if(nrow(top10_nodes)>0)fwrite(top10_nodes,f))
}

shinyApp(ui=ui, server=server)
