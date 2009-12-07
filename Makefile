#!/usr/bin/make

ifneq ($(CONFIG),)
include Makefile.config.$(CONFIG)
else
include Makefile.config
endif

.DEFAULT_GOAL= all

FIRELIBS =  $(LIBDIR)/libRetriever.a $(LIBDIR)/libClustering.a  $(LIBDIR)/libDistanceFunctions.a $(LIBDIR)/libFeatureExtractors.a $(LIBDIR)/libFeatures.a  $(LIBDIR)/libImage.a $(LIBDIR)/libCore.a 

# Core ------------------------------------------------------------
LIBCORE_SOURCES = Core/diag.cpp Core/gzstream.cpp Core/hungarian.cpp Core/jflib.cpp Core/Lapack.cpp Core/lda.cpp Core/pca.cpp Core/runprogram.cpp Core/ScopeTimer.cpp Core/svd.cpp Core/net.cpp Core/supportvectormachine.cpp Core/stringparser.cpp
LIBCORE_OBJECTS := $(patsubst %.o,$(OBJDIR)/%.o,$(LIBCORE_SOURCES:.cpp=.o))
$(LIBDIR)/libCore.a: $(LIBCORE_OBJECTS)

# Distances--------------------------------------------------------
LIBDISTANCES_SOURCES =    Retriever/getscoring.cpp Retriever/maxentscoring.cpp Retriever/maxentscoringfirstandsecondorder.cpp Retriever/maxentscoringsecondorder.cpp Retriever/distancemaker.cpp Retriever/distancemaker.cpp DistanceFunctions/dist_distfile.cpp DistanceFunctions/dist_bm25.cpp DistanceFunctions/dist_globallocalfeaturedistance.cpp DistanceFunctions/dist_idm.cpp DistanceFunctions/dist_lfhungarian.cpp DistanceFunctions/dist_lfsigemd.cpp DistanceFunctions/dist_metafeature.cpp DistanceFunctions/dist_mpeg7.cpp DistanceFunctions/dist_rast.cpp DistanceFunctions/dist_smart2.cpp DistanceFunctions/dist_textfeature.cpp DistanceFunctions/dist_tfidf.cpp DistanceFunctions/emd.cpp DistanceFunctions/dist_weightedl1.cpp	
LIBDISTANCES_OBJECTS := $(patsubst %.o,$(OBJDIR)/%.o,$(LIBDISTANCES_SOURCES:.cpp=.o))
$(LIBDIR)/libDistanceFunctions.a: $(LIBDISTANCES_OBJECTS)

# Retriever -------------------------------------------------------
LIBRETRIEVER_SOURCES = Retriever/database.cpp     Retriever/featureloader.cpp  Retriever/imagecomparator.cpp  Retriever/largebinaryfeaturefile.cpp Retriever/largefeaturefile.cpp  Retriever/retriever.cpp Retriever/server.cpp Retriever/querycombiner.cpp Retriever/reranker.cpp
LIBRETRIEVER_OBJECTS := $(patsubst %.o,$(OBJDIR)/%.o,$(LIBRETRIEVER_SOURCES:.cpp=.o))
$(LIBDIR)/libRetriever.a: $(LIBRETRIEVER_OBJECTS)

RETRIEVER_SOURCES = Retriever/fire.cpp
RETRIEVER_OBJECTS :=  $(patsubst %.o,$(OBJDIR)/%.o,$(RETRIEVER_SOURCES:.cpp=.o))
RETRIEVER_PROGRAMS := $(patsubst $(OBJDIR)/Retriever/%.o,$(BINDIR)/%,$(RETRIEVER_OBJECTS))
$(BINDIR)/fire: $(OBJDIR)/Retriever/fire.o $(FIRELIBS)

# Image -------------------------------------------------------
LIBIMAGE_SOURCES = Image/colorhsv.cpp Image/imagelib.cpp Image/interpolatingimage.cpp
LIBIMAGE_OBJECTS := $(patsubst %.o,$(OBJDIR)/%.o,$(LIBIMAGE_SOURCES:.cpp=.o))
$(LIBDIR)/libImage.a: $(LIBIMAGE_OBJECTS)

# FeatureExtractors -------------------------------------------------------
LIBFEATEX_SOURCES = FeatureExtractors/createsparsehisto.cpp FeatureExtractors/differenceofgaussian.cpp FeatureExtractors/gabor.cpp FeatureExtractors/globalfeatureextraction.cpp FeatureExtractors/invariantfeaturehistogram.cpp FeatureExtractors/kernelfunctionmaker.cpp FeatureExtractors/localfeatureextractor.cpp FeatureExtractors/relationalfeaturehistogram.cpp FeatureExtractors/salientpoints.cpp FeatureExtractors/extracttemplate.cpp FeatureExtractors/sift.cpp FeatureExtractors/tamurafeature.cpp FeatureExtractors/wavelet.cpp
LIBFEATEX_OBJECTS := $(patsubst %.o,$(OBJDIR)/%.o,$(LIBFEATEX_SOURCES:.cpp=.o))
$(LIBDIR)/libFeatureExtractors.a: $(LIBFEATEX_OBJECTS)
FEATEX_SOURCES = FeatureExtractors/extractlocalfeatures.cpp FeatureExtractors/extractrelationalfeaturehistogram.cpp FeatureExtractors/extractrelationaltexturefeaturepairs.cpp FeatureExtractors/extractsift.cpp FeatureExtractors/extractsparsecolorhistogram.cpp FeatureExtractors/extractsparsegaborhistogram.cpp FeatureExtractors/extractsparsepatchhistogram.cpp FeatureExtractors/extractsparsesifthistogram.cpp FeatureExtractors/extractsparsesurfhistogram.cpp FeatureExtractors/extracttamuratexturefeature.cpp FeatureExtractors/extracttamuratexturefeaturepairs.cpp FeatureExtractors/calcsalientpoints.cpp FeatureExtractors/extractaspectratio.cpp FeatureExtractors/extractcolorhistogram.cpp FeatureExtractors/extractgaborfeaturevector.cpp FeatureExtractors/extractglobaltexturefeature.cpp FeatureExtractors/extractinvariantfeaturehistogram.cpp FeatureExtractors/extractinvarianttexturefeaturepairs.cpp FeatureExtractors/colororgray.cpp  
FEATEX_OBJECTS := $(patsubst %.o,$(OBJDIR)/%.o,$(FEATEX_SOURCES:.cpp=.o))
FEATEX_PROGRAMS := $(patsubst FeatureExtractors/%.o,$(BINDIR)/%,$(FEATEX_SOURCES:.cpp=.o))
$(BINDIR)/extractcolorhistogram:$(OBJDIR)/FeatureExtractors/extractcolorhistogram.o $(FIRELIBS)
$(BINDIR)/extractlocalfeatures:$(OBJDIR)/FeatureExtractors/extractlocalfeatures.o $(FIRELIBS)
$(BINDIR)/extractlf:$(OBJDIR)/FeatureExtractors/extractlf.o $(FIRELIBS)
$(BINDIR)/extractrelationalfeaturehistogram:$(OBJDIR)/FeatureExtractors/extractrelationalfeaturehistogram.o $(FIRELIBS)
$(BINDIR)/extractsift:$(OBJDIR)/FeatureExtractors/extractsift.o $(FIRELIBS)
$(BINDIR)/extractlf: $(OBJDIR)/FeatureExtractors/extractlf.o $(FIRELIBS)
$(BINDIR)/extractlocalfeatures: $(OBJDIR)/FeatureExtractors/extractlocalfeatures.o $(FIRELIBS)
$(BINDIR)/extractrelationalfeaturehistogram: $(OBJDIR)/FeatureExtractors/extractrelationalfeaturehistogram.o $(FIRELIBS)
$(BINDIR)/extractrelationaltexturefeaturepairs: $(OBJDIR)/FeatureExtractors/extractrelationaltexturefeaturepairs.o $(FIRELIBS)
$(BINDIR)/extractsift: $(OBJDIR)/FeatureExtractors/extractsift.o $(FIRELIBS)
$(BINDIR)/extractsparsecolorhistogram: $(OBJDIR)/FeatureExtractors/extractsparsecolorhistogram.o $(FIRELIBS)
$(BINDIR)/extractsparsegaborhistogram: $(OBJDIR)/FeatureExtractors/extractsparsegaborhistogram.o $(FIRELIBS)
$(BINDIR)/extractsparsepatchhistogram: $(OBJDIR)/FeatureExtractors/extractsparsepatchhistogram.o $(FIRELIBS)
$(BINDIR)/extractsparsesifthistogram: $(OBJDIR)/FeatureExtractors/extractsparsesifthistogram.o $(FIRELIBS)
$(BINDIR)/extractsparsesurfhistogram: $(OBJDIR)/FeatureExtractors/extractsparsesurfhistogram.o $(FIRELIBS)
$(BINDIR)/extracttamuratexturefeature: $(OBJDIR)/FeatureExtractors/extracttamuratexturefeature.o $(FIRELIBS)
$(BINDIR)/extracttamuratexturefeaturepairs: $(OBJDIR)/FeatureExtractors/extracttamuratexturefeaturepairs.o $(FIRELIBS)
$(BINDIR)/calcsalientpoints: $(OBJDIR)/FeatureExtractors/calcsalientpoints.o $(FIRELIBS)
$(BINDIR)/extractaspectratio: $(OBJDIR)/FeatureExtractors/extractaspectratio.o $(FIRELIBS)
$(BINDIR)/extractcolorhistogram: $(OBJDIR)/FeatureExtractors/extractcolorhistogram.o $(FIRELIBS)
$(BINDIR)/extractgaborfeaturevector: $(OBJDIR)/FeatureExtractors/extractgaborfeaturevector.o $(FIRELIBS)
$(BINDIR)/extractglobaltexturefeature: $(OBJDIR)/FeatureExtractors/extractglobaltexturefeature.o $(FIRELIBS)
$(BINDIR)/extractinvariantfeaturehistogram: $(OBJDIR)/FeatureExtractors/extractinvariantfeaturehistogram.o $(FIRELIBS)
$(BINDIR)/colororgray: $(OBJDIR)/FeatureExtractors/colororgray.o $(FIRELIBS)
$(BINDIR)/extractinvarianttexturefeaturepairs: $(OBJDIR)/FeatureExtractors/extractinvarianttexturefeaturepairs.o $(FIRELIBS)
$(BINDIR)/lftolfsignature: $(OBJDIR)/FeatureExtractors/lftolfsignature.o $(FIRELIBS)

# Features -------------------------------------------------------
LIBFEATURES_SOURCES = Features/basefeature.cpp Features/histogramfeature.cpp Features/imagefeature.cpp Features/lfsignaturefeature.cpp Features/localfeatures.cpp Features/mpeg7feature.cpp Features/sparsehistogramfeature.cpp Features/vectorfeature.cpp Features/lfposclusteridfeature.cpp Features/imagecontainer.cpp
LIBFEATURES_OBJECTS := $(patsubst %.o,$(OBJDIR)/%.o,$(LIBFEATURES_SOURCES:.cpp=.o))
$(LIBDIR)/libFeatures.a: $(LIBFEATURES_OBJECTS)


# Clustering ---------------------------------------------------------
LIBCLUSTERING_SOURCES = Clustering/clusterlocalfeatures.cpp Clustering/dbscan.cpp Clustering/em.cpp Clustering/gmd.cpp  Clustering/positionclusterer.cpp
LIBCLUSTERING_OBJECTS := $(patsubst %.o,$(OBJDIR)/%.o,$(LIBCLUSTERING_SOURCES:.cpp=.o))
$(LIBDIR)/libClustering.a: $(LIBCLUSTERING_OBJECTS)
CLUSTERING_SOURCES = Clustering/emclustercenter2vectorfeature.cpp Clustering/imageclusterer.cpp Clustering/jfclustering.cpp Clustering/lfclustering.cpp  Clustering/visualizeclusterpositions.cpp
CLUSTERING_OBJECTS := $(patsubst %.o,$(OBJDIR)/%.o,$(CLUSTERING_SOURCES:.cpp=.o))
CLUSTERING_PROGRAMS := $(patsubst Clustering/%.o,$(BINDIR)/%,$(CLUSTERING_SOURCES:.cpp=.o))
$(BINDIR)/emclustercenter2vectorfeature: $(OBJDIR)/Clustering/emclustercenter2vectorfeature.o $(FIRELIBS)
$(BINDIR)/imageclusterer: $(OBJDIR)/Clustering/imageclusterer.o $(FIRELIBS)
$(BINDIR)/jfclustering: $(OBJDIR)/Clustering/jfclustering.o $(FIRELIBS)
$(BINDIR)/lfclustering: $(OBJDIR)/Clustering/lfclustering.o $(FIRELIBS)
$(BINDIR)/visualizeclusterpositions: $(OBJDIR)/Clustering/visualizeclusterpositions.o $(FIRELIBS)

#Tools ----------------------------------------------------------------
TOOLS_SOURCES = Tools/db2jf.cpp Tools/db2lbff.cpp Tools/db2lff.cpp Tools/gaborcreatejf.cpp Tools/jf2arff.cpp Tools/lfcreatejf.cpp Tools/mergejf.cpp Tools/pixelgmd.cpp Tools/randomfilelistcreator.cpp Tools/sobelcreatejf.cpp Tools/sparsehisto2histo.cpp
TOOLS_OBJECTS := $(patsubst %.o,$(OBJDIR)/%.o,$(TOOLS_SOURCES:.cpp=.o))
TOOLS_PROGRAMS := $(patsubst Tools/%.o,$(BINDIR)/%,$(TOOLS_SOURCES:.cpp=.o))
$(BINDIR)/db2jf: $(OBJDIR)/Tools/db2jf.o $(FIRELIBS)
$(BINDIR)/db2lbff: $(OBJDIR)/Tools/db2lbff.o $(FIRELIBS)
$(BINDIR)/db2lff: $(OBJDIR)/Tools/db2lff.o $(FIRELIBS)
$(BINDIR)/gaborcreatejf: $(OBJDIR)/Tools/gaborcreatejf.o $(FIRELIBS)
$(BINDIR)/jf2arff: $(OBJDIR)/Tools/jf2arff.o $(FIRELIBS)
$(BINDIR)/lfcreatejf: $(OBJDIR)/Tools/lfcreatejf.o $(FIRELIBS)
$(BINDIR)/mergejf: $(OBJDIR)/Tools/mergejf.o $(FIRELIBS)
$(BINDIR)/pixelgmd: $(OBJDIR)/Tools/pixelgmd.o $(FIRELIBS)
$(BINDIR)/randomfilelistcreator: $(OBJDIR)/Tools/randomfilelistcreator.o $(FIRELIBS)
$(BINDIR)/sobelcreatejf: $(OBJDIR)/Tools/sobelcreatejf.o $(FIRELIBS)
$(BINDIR)/sparsehisto2histo: $(OBJDIR)/Tools/sparsehisto2histo.o $(FIRELIBS)
$(BINDIR)/lf2png: $(OBJDIR)/Tools/lf2png.o $(FIRELIBS)
$(BINDIR)/findduplicates: $(OBJDIR)/Tools/findduplicates.o $(FIRELIBS)

# Misc ----------------------------------------------------------------
MISC_SOURCES = Misc/collage.cpp Misc/dbpca.cpp Misc/facefeatureprocessor.cpp Misc/featurescomparator.cpp Misc/eigenfacer.cpp Misc/histogramnormalization.cpp Misc/mosaic.cpp  Misc/pcavectortoimage.cpp Misc/testscaleinvariantfeatures.cpp Misc/testsparsehistogramfeature.cpp Misc/visualizelocalfeatures.cpp 
MISC_OBJECTS := $(patsubst %.o,$(OBJDIR)/%.o,$(MISC_SOURCES:.cpp=.o))
MISC_PROGRAMS := $(patsubst Misc/%.o,$(BINDIR)/%,$(MISC_SOURCES:.cpp=.o))
$(BINDIR)/collage: $(OBJDIR)/Misc/collage.o $(FIRELIBS)
$(BINDIR)/create_classes: $(OBJDIR)/Misc/create_classes.o $(FIRELIBS)
$(BINDIR)/dbpca: $(OBJDIR)/Misc/dbpca.o $(FIRELIBS)
$(BINDIR)/eigenfacer: $(OBJDIR)/Misc/eigenfacer.o $(FIRELIBS)
$(BINDIR)/facefeatureprocessor: $(OBJDIR)/Misc/facefeatureprocessor.o $(FIRELIBS)
$(BINDIR)/featurescomparator: $(OBJDIR)/Misc/featurescomparator.o $(FIRELIBS)
$(BINDIR)/histogramnormalization: $(OBJDIR)/Misc/histogramnormalization.o $(FIRELIBS)
$(BINDIR)/mosaic: $(OBJDIR)/Misc/mosaic.o $(FIRELIBS)
$(BINDIR)/pcavectortoimage: $(OBJDIR)/Misc/pcavectortoimage.o $(FIRELIBS)
$(BINDIR)/testscaleinvariantfeatures: $(OBJDIR)/Misc/testscaleinvariantfeatures.o $(FIRELIBS)
$(BINDIR)/testsparsehistogramfeature: $(OBJDIR)/Misc/testsparsehistogramfeature.o $(FIRELIBS)
$(BINDIR)/vis-rast-matching: $(OBJDIR)/Misc/vis-rast-matching.o $(FIRELIBS)
$(BINDIR)/visualizelocalfeatures: $(OBJDIR)/Misc/visualizelocalfeatures.o $(FIRELIBS)

all: $(BINDIR) $(OBJDIR) $(LIBDIR) $(RETRIEVER_PROGRAMS) $(FEATEX_PROGRAMS) $(CLASSIFIERS_PROGRAMS) $(CLUSTERING_PROGRAMS) $(TEXTIR_PROGRAMS) $(MISC_PROGRAMS) $(TOOLS_PROGRAMS)

info:
	@echo $(FEATEX_PROGRAMS)
	@echo $(RETRIEVER_PROGRAMS)
	@echo $(DEPFILES)
	@echo $(ARCH)

$(BINDIR)/%: 
	mkdir -p $(BINDIR)
	$(LD) $(filter %.o %.a ,$^) -o $@ $(LDFLAGS)

DEPFILES = $(LIBCORE_OBJECTS:.o=.d) $(LIBDISTANCES_OBJECTS:.o=.d) $(LIBRETRIEVER_OBJECTS:.o=.d) $(RETRIEVER_OBJECTS:.o=.d) $(LIBIMAGE_OBJECTS:.o=.d) $(LIBFEATEX_OBJECTS:.o=.d) $(FEATEX_OBJECTS:.o=.d) $(LIBFEATURES_OBJECTS:.o=.d) $(CLASSIFIERS_OBJECTS:.o=.d) $(LIBCLUSTERING_OBJECTS:.o=.d) $(CLUSTERING_OBJECTS:.o=.d) $(MISC_OBJECTS:.o=.d) $(TOOLS_OBJECTS:.o=.d) 
-include $(DEPFILES)

include Makefile.rules
