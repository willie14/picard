package picard.arrays.illumina;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.samtools.reference.ReferenceSequenceFileFactory;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalList;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.ProgressLogger;
import htsjdk.tribble.annotation.Strand;
import htsjdk.variant.utils.SAMSequenceDictionaryExtractor;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;
import org.apache.commons.lang.StringUtils;
import org.apache.commons.lang3.ArrayUtils;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import picard.PicardException;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.StandardOptionDefinitions;
import picard.vcf.ByIntervalListVariantContextIterator;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileOutputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

/**
 * Create an Extended Illumina Manifest by performing a liftover to Build 37.
 */
@CommandLineProgramProperties(
        summary = "Create an Extended Illumina Manifest",
        oneLineSummary = "Create an Extended Illumina Manifest by performing a liftover to Build 37",
        programGroup = picard.cmdline.programgroups.GenotypingArraysProgramGroup.class
)
public class CreateExtendedIlluminaManifest extends CommandLineProgram {

    @Argument(shortName = StandardOptionDefinitions.INPUT_SHORT_NAME, doc = "This is the text version of the Illumina .bpm file")
    public File INPUT;

    @Argument(shortName = "BPM", doc = "The Illumina Bead Pool Manifest (.bpm) file")
    public File BEAD_POOL_MANIFEST_FILE;

    @Argument(shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc = "The name of the extended manifest to be written.")
    public File OUTPUT;

    @Argument(shortName = "BAF", doc = "The name of the the 'bad assays file'. This is a subset version of the extended manifest, " +
            "containing only unmappable assays")
    public File BAD_ASSAYS_FILE;

    @Argument(shortName = "RF", doc = "The name of the the report file")
    public File REPORT_FILE;

    @Argument(shortName = "CF", doc = "The Standard (Hapmap-trained) cluster file (.egt) from Illumina. " +
            "If there are duplicate assays at a site, this is used to decide which is the 'best' (non-filtered in generated VCFs) " +
            "by choosing the assay with the best GenTrain scores)")
    public File CLUSTER_FILE;

    @Argument(shortName = "DBSNP", doc = "Reference dbSNP file in VCF format.", optional = true)
    public File DBSNP_FILE;

    @Argument(shortName = "TB", doc = "The target build.")
    public String TARGET_BUILD;

    @Argument(shortName = "TR", doc = "The target build's reference file.")
    public File TARGET_REFERENCE_FILE;

    @Argument(shortName = "SB", doc = "A supported build. The order of the input must match the order for SUPPORTED_REFERENCE_FILE and SUPPORTED_CHAIN_FILE.", optional = true)
    public List<String> SUPPORTED_BUILD;

    @Argument(shortName = "SR", doc = "A reference file for a supported build. Must provide a supported chain file to convert from supported -> target.", optional = true)
    public List<File> SUPPORTED_REFERENCE_FILE;

    @Argument(shortName = "SC", doc = "A chain file that maps from a supported build -> target build. Must provide a corresponding supported reference file.", optional = true)
    public List<File> SUPPORTED_CHAIN_FILE;

    private static final Log log = Log.getInstance(CreateExtendedIlluminaManifest.class);

    public static final String VERSION = "1.6";

    @Override
    protected int doWork() {

        try {
            // Load the sequence dictionary from the Target Reference file
            final SAMSequenceDictionary sequenceDictionary = SAMSequenceDictionaryExtractor.extractDictionary(TARGET_REFERENCE_FILE);

            ProgressLogger logger = new ProgressLogger(log, 10000);
            final Map<String, ReferenceSequenceFile> referenceSequenceMap = new HashMap<>();
            final Map<String, File> chainFilesMap = new HashMap<>();

            referenceSequenceMap.put(TARGET_BUILD, ReferenceSequenceFileFactory.getReferenceSequenceFile(TARGET_REFERENCE_FILE));

            for (int i = 0; i < SUPPORTED_BUILD.size(); i++) {
                referenceSequenceMap.put(SUPPORTED_BUILD.get(i), ReferenceSequenceFileFactory.getReferenceSequenceFile(SUPPORTED_REFERENCE_FILE.get(i)));
                chainFilesMap.put(SUPPORTED_BUILD.get(i), SUPPORTED_CHAIN_FILE.get(i));
            }

            // Open the Original Illumina Manifest
            final IlluminaManifest manifestFile = new IlluminaManifest(INPUT);

            // TODO - Warn - don't allow overwrite
            IOUtil.assertFileIsWritable(OUTPUT);
            IOUtil.assertFileIsWritable(BAD_ASSAYS_FILE);
            IOUtil.assertFileIsWritable(REPORT_FILE);

            Map<String, List<IlluminaManifestRecord>> coordinateMap = new HashMap<>();

            IntervalList manifestSnpIntervals = new IntervalList(sequenceDictionary);
            IntervalList manifestIndelIntervals = new IntervalList(sequenceDictionary);

            log.info("Loading the bpm file");
            final IlluminaBPMFile illuminaBPMFile;
            try {
                illuminaBPMFile = new IlluminaBPMFile(BEAD_POOL_MANIFEST_FILE);
            } catch (IOException e) {
                throw new PicardException("Error reading bpm file '" + BEAD_POOL_MANIFEST_FILE.getAbsolutePath() + "'", e);
            }
            IlluminaBPMLocusEntry[] illuminaBPMLocusEntries = illuminaBPMFile.getLocusEntries();

            // Load the cluster file to get the GenTrain scores
            // load the egt first, and create a map of ilmnid to gentrain score.  Save that and use it for deduplicating.
            log.info("Loading the egt file");
            final InfiniumEGTFile infiniumEGTFile;
            try {
                infiniumEGTFile = new InfiniumEGTFile(CLUSTER_FILE);
            } catch (IOException e) {
                throw new PicardException("Error reading cluster file '" + CLUSTER_FILE.getAbsolutePath() + "'", e);
            }

            Build37ExtendedIlluminaManifestRecordCreator creator = new Build37ExtendedIlluminaManifestRecordCreator(referenceSequenceMap, chainFilesMap);

            // first iteration through the manifest to find all dupes
            log.info("Phase 1.  First Pass through the manifest.  Build coordinate map for dupe flagging and make SNP and indel-specific interval lists for parsing dbSnp");
            final Iterator<IlluminaManifestRecord> firstPassIterator = manifestFile.iterator();

            ManifestStatistics manifestStatistics = new ManifestStatistics();

            //grab the dictionary from the VCF and use it in the IntervalList
            final SAMFileHeader samFileHeader = new SAMFileHeader();
            samFileHeader.setSequenceDictionary(sequenceDictionary);

            int locusIndex = 0;
            while (firstPassIterator.hasNext()) {
                logger.record("0", 0);
                if (locusIndex > illuminaBPMLocusEntries.length) {
                    throw new PicardException("Differing number of entries between bpm and manifest file");
                }
                IlluminaBPMLocusEntry locusEntry = illuminaBPMLocusEntries[locusIndex++];
                final IlluminaManifestRecord record = firstPassIterator.next();

                creator.validateLocusEntryAndCreateExtendedRecord(locusEntry, record);

                // Create an ExtendedIlluminaManifestRecord here so that we can get the (potentially lifted over) coordinates
                final Build37ExtendedIlluminaManifestRecord rec = new Build37ExtendedIlluminaManifestRecord(record,
                        referenceSequenceMap, chainFilesMap, false, null);
                manifestStatistics.updateStatistics(rec);

                // A DUP is only a DUP if it's at the same location AND has the same alleles...
                String key = rec.getB37Chr() + ":" + rec.getB37Pos() + "." + rec.getAlleleA().toString() + "." + rec.getAlleleB();
                if (coordinateMap.containsKey(key)) {
                    coordinateMap.get(key).add(record);
                } else {
                    List<IlluminaManifestRecord> newList = new ArrayList<>();
                    newList.add(record);
                    coordinateMap.put(key, newList);
                }

                if (!rec.isBad()) {
                    final int length = Integer.max(rec.getAlleleA().length(), rec.getAlleleB().length());
                    Interval interval = new Interval(rec.getB37Chr(), rec.getB37Pos(), rec.getB37Pos() + length);
                    if (rec.isSnp()) {
                        manifestSnpIntervals.add(interval);
                    } else {
                        manifestIndelIntervals.add(interval);
                    }
                }
            }

            // Generate a sorted set of the variants in the Illumina Manifest so that we can check them
            // Against the (sorted) dbSnp Vcf.
            manifestSnpIntervals = manifestSnpIntervals.sorted();
            manifestIndelIntervals = manifestIndelIntervals.sorted();

            log.info("Phase 2.  Parse dbSnpVCF and build SNP and indel-specific locus to rsId maps");
            Map<String, String> snpLocusToRsId = new HashMap<>();
            Map<String, String> indelLocusToRsId = new HashMap<>();
            if (DBSNP_FILE != null) {
                // Because dbSnp can contain both SNPs and indels which may be at the same locus,
                // We do two passes through dbSnpVcf to build separate maps.
                log.info("SNP-specific");
                snpLocusToRsId = generateLocusToRsidMap(DBSNP_FILE, manifestSnpIntervals);
                log.info("indel-specific");
                indelLocusToRsId = generateLocusToRsidMap(DBSNP_FILE, manifestIndelIntervals);
            }

            // filter out all unique coordinates
            Map<String, List<IlluminaManifestRecord>> dupeMap = coordinateMap.entrySet().stream()
                    .filter(map -> map.getValue().size() > 1)
                    .collect(Collectors.toMap(Map.Entry::getKey, Map.Entry::getValue));
            coordinateMap.clear();

            // evaluate each coordinate assay and remove the assay with the best GenTrain score (all remaining are dupes)
            dupeMap.entrySet().forEach(entry ->
                    entry.getValue().remove(entry.getValue().stream().max(Comparator.comparingDouble(assay ->
                            infiniumEGTFile.totalScore[infiniumEGTFile.rsNameToIndex.get(assay.getName())])).get()));

            // we really only need the list of indices for the dupes
            List<Integer> dupeIndices = dupeMap.entrySet().stream()
                    .flatMapToInt(entry ->
                            entry.getValue().stream()
                                    .mapToInt(IlluminaManifestRecord::getIndex))
                    .boxed().collect(Collectors.toList());
            dupeMap.clear();

            final IlluminaManifest secondPassManifestFile = new IlluminaManifest(INPUT);
            final Iterator<IlluminaManifestRecord> secondPassIterator = secondPassManifestFile.iterator();

            final BufferedWriter out = new BufferedWriter(new FileWriter(OUTPUT, false));
            writeExtendedIlluminaManifestHeaders(manifestFile, out);

            //second iteration to write all records after dupe evaluation
            log.info("Phase 3.  Generate the Extended Illumina Manifest");
            logger = new ProgressLogger(log, 10000);
            List<Build37ExtendedIlluminaManifestRecord> badRecords = new ArrayList<>();
            while (secondPassIterator.hasNext()) {
                logger.record("0", 0);
                final IlluminaManifestRecord record = secondPassIterator.next();
                final String locus = record.getChr() + "." + record.getPosition();
                String rsId;
                if (record.isSnp()) {
                    rsId = snpLocusToRsId.get(locus);
                } else {
                    rsId = indelLocusToRsId.get(locus);
                }
                final Build37ExtendedIlluminaManifestRecord rec = new Build37ExtendedIlluminaManifestRecord(record,
                        referenceSequenceMap, chainFilesMap, dupeIndices.contains(record.getIndex()), rsId);
                if (rec.isBad()) {
                    badRecords.add(rec);
                }
                out.write(rec.getLine());
                out.newLine();
            }

            out.flush();
            out.close();

            writeBadAssaysFile(BAD_ASSAYS_FILE, badRecords);
            manifestStatistics.logStatistics(REPORT_FILE);
        } catch (IOException e) {
            throw new PicardException(e.getMessage(), e);
        }

        return 0;
    }

    private void writeBadAssaysFile(File badAssaysFile, List<Build37ExtendedIlluminaManifestRecord> badRecords) throws IOException {
        BufferedWriter badAssaysFileWriter;
        badAssaysFileWriter = new BufferedWriter(new FileWriter(badAssaysFile, false));
        badAssaysFileWriter.write("## The following assays were marked by CreateExtendedIlluminaManifest as Unparseable (input file: " + INPUT.getAbsolutePath() + ")");
        badAssaysFileWriter.newLine();
        badAssaysFileWriter.write("#IlmnId,Name,GenomeBuild,Chr,MapInfo,FailureFlag");
        badAssaysFileWriter.newLine();

        for (Build37ExtendedIlluminaManifestRecord record : badRecords) {
            final List<String> badRecord = java.util.Arrays.asList(record.getIlmnId(), record.getName(), record.getGenomeBuild(), record.getChr(), "" + record.getPosition(), record.getFlag().toString());
            badAssaysFileWriter.write(StringUtils.join(badRecord, ","));
            badAssaysFileWriter.newLine();
        }
        badAssaysFileWriter.flush();
        badAssaysFileWriter.close();
    }

    private static class ManifestStatistics {
        int numAssays;
        int numSnps;
        int numIndels;
        int numAssaysFlagged;
        int numSnpsFlagged;
        int numSnpProbeSequenceMismatch;
        int numAmbiguousSnpsOnPosStrand;
        int numAmbiguousSnpsOnNegStrand;

        int numIndelsFlagged;
        int numIndelProbeSequenceMismatch;
        int numIndelProbeSequenceStrandInvalid;
        int numIndelSourceSequenceMismatch;
        int numIndelSourceSequenceStrandInvalid;
        int numIndelSourceSequenceInvalid;
        int numIndelsNotFound;
        int numIndelConfict;

        int numOnBuild37;
        int numOnBuild36;
        int numOnOtherBuild;
        int numLiftoverFailed;

        void updateStatistics(Build37ExtendedIlluminaManifestRecord rec) {
            numAssays++;
            if (rec.isSnp()) {
                numSnps++;
            } else {
                numIndels++;
            }
            if (rec.getMajorGenomeBuild().equals(Build37ExtendedIlluminaManifestRecord.BUILD_37)) {
                numOnBuild37++;
            } else if (rec.getMajorGenomeBuild().equals(Build37ExtendedIlluminaManifestRecord.BUILD_36)) {
                numOnBuild36++;
            } else {
                numOnOtherBuild++;
            }
            if (rec.getFlag().equals(Build37ExtendedIlluminaManifestRecord.Flag.LIFTOVER_FAILED)) {
                numLiftoverFailed++;
            }
            if (!rec.isBad()) {
                if (rec.isAmbiguous()) {
                    if (rec.getCalculatedStrand() == Strand.NEGATIVE) {
                        numAmbiguousSnpsOnNegStrand++;
                    }
                    if (rec.getCalculatedStrand() == Strand.POSITIVE) {
                        numAmbiguousSnpsOnPosStrand++;
                    }
                }
            } else {
                numAssaysFlagged++;
                if (rec.isSnp()) {
                    numSnpsFlagged++;
                } else {
                    numIndelsFlagged++;
                }
                if (rec.isIndel()) {
                    switch (rec.getFlag()) {
                        case PROBE_SEQUENCE_MISMATCH:
                            numIndelProbeSequenceMismatch++;
                            break;
                        case PROBE_SEQUENCE_STRAND_INVALID:
                            numIndelProbeSequenceStrandInvalid++;
                            break;
                        case SOURCE_SEQUENCE_MISMATCH:
                            numIndelSourceSequenceMismatch++;
                            break;
                        case SOURCE_SEQUENCE_INVALID:
                            numIndelSourceSequenceInvalid++;
                            break;
                        case SOURCE_SEQUENCE_STRAND_INVALID:
                            numIndelSourceSequenceStrandInvalid++;
                            break;
                        case INDEL_NOT_FOUND:
                            numIndelsNotFound++;
                            break;
                        case INDEL_CONFLICT:
                            numIndelConfict++;
                            break;
                    }
                }
                else {
                    if (rec.getFlag() == Build37ExtendedIlluminaManifestRecord.Flag.PROBE_SEQUENCE_MISMATCH) {
                        numSnpProbeSequenceMismatch++;
                    }
                }
            }
        }

        void logStatistics(File output) throws IOException {
            try (BufferedWriter writer = new BufferedWriter(new OutputStreamWriter(
                    new FileOutputStream(output), "utf-8"))) {
                writer.write("Number of assays written: " + numAssays);
                writer.newLine();
                writer.write("Number of assays flagged: " + numAssaysFlagged);
                writer.newLine();
                writer.write("Number of SNPs: " + numSnps);
                writer.newLine();
                writer.write("Number of SNPs flagged: " + numSnpsFlagged);
                writer.newLine();
                writer.write("Number of SNPs flagged for sequence mismatch: " + numSnpProbeSequenceMismatch);
                writer.newLine();
                writer.write("Number of ambiguous SNPs on Positive Strand: " + numAmbiguousSnpsOnPosStrand);
                writer.newLine();
                writer.write("Number of ambiguous SNPs on Negative Strand: " + numAmbiguousSnpsOnNegStrand);
                writer.newLine();

                writer.write("Number of Indels: " + numIndels);
                writer.newLine();
                writer.write("Number of Indels flagged: " + numIndelsFlagged);
                writer.newLine();
                writer.write("Number of Indels flagged for probe sequence mismatch: " + numIndelProbeSequenceMismatch);
                writer.newLine();
                writer.write("Number of Indels flagged for probe sequence strand invalid: " + numIndelProbeSequenceStrandInvalid);
                writer.newLine();
                writer.write("Number of Indels flagged for source sequence mismatch: " + numIndelSourceSequenceMismatch);
                writer.newLine();
                writer.write("Number of Indels flagged for source sequence invalid: " + numIndelSourceSequenceInvalid);
                writer.newLine();
                writer.write("Number of Indels flagged for source sequence strand invalid: " + numIndelSourceSequenceStrandInvalid);
                writer.newLine();
                writer.write("Number of Indels not found: " + numIndelsNotFound);
                writer.newLine();
                writer.write("Number of Indels flagged for conflict: " + numIndelConfict);
                writer.newLine();

                writer.write("Number of assays on Build37: " + numOnBuild37);
                writer.newLine();
                writer.write("Number of assays on Build36: " + numOnBuild36);
                writer.newLine();
                writer.write("Number of assays on Other Build: " + numOnOtherBuild);
                writer.newLine();
                writer.write("Number of assays failing liftover: " + numLiftoverFailed);
                writer.newLine();
            }
        }
    }

    @Override
    protected String[] customCommandLineValidation() {
        IOUtil.assertFileIsReadable(INPUT);
        IOUtil.assertFileIsReadable(CLUSTER_FILE);
        if (DBSNP_FILE != null) {
            IOUtil.assertFileIsReadable(DBSNP_FILE);
        }

        IOUtil.assertFileIsReadable(TARGET_REFERENCE_FILE);
        for (File f : SUPPORTED_REFERENCE_FILE) IOUtil.assertFileIsReadable(f);
        for (File f : SUPPORTED_CHAIN_FILE) IOUtil.assertFileIsReadable(f);

        final List<String> errors = new ArrayList<>();

        if (SUPPORTED_BUILD.size() != SUPPORTED_REFERENCE_FILE.size()) {
            errors.add("The number of supported builds does not match the number of supported reference files");
        }

        if (SUPPORTED_BUILD.size() != SUPPORTED_CHAIN_FILE.size()) {
            errors.add("The number of supported builds does not match the number of supported chain files");
        }

        return (errors.size() > 0)
                ? errors.toArray(new String[errors.size()])
                : null;
    }

    /**
     * Generates a mapping of locus (contig.posn) in the manifest file to rsId.
     * Uses the passed interval list to selectively parse dbSnpVcf.
     * Returns a map of locus to rsId
     * @param dbSnpFile the dbSnp file to parse.
     * @param intervals interval list for which intervals to parse out of the dbSnp file
     * @return mapping of locus in the manifest file (contig.posn) to rsId
     */
    Map<String, String> generateLocusToRsidMap(File dbSnpFile, IntervalList intervals) {
        ProgressLogger logger = new ProgressLogger(log, 10000);

        Map<String, String> manifestLocusToRsId = new HashMap<>();

        final VCFFileReader dbSnpReader = new VCFFileReader(dbSnpFile, true);
        final Iterator<VariantContext> dbSnpIterator = new ByIntervalListVariantContextIterator(dbSnpReader, intervals);
        while (dbSnpIterator.hasNext()) {
            VariantContext variantContext = dbSnpIterator.next();
            logger.record(variantContext.getContig(), variantContext.getStart());

            for (int posn = variantContext.getStart(); posn <= variantContext.getEnd(); posn++) {
                final String locus = variantContext.getContig() + "." + posn;
                manifestLocusToRsId.put(locus, variantContext.getID());
            }
        }

        return manifestLocusToRsId;
    }



    void writeExtendedIlluminaManifestHeaders(final IlluminaManifest manifest, final BufferedWriter output) throws IOException {
        int numColumns = -1;
        List<String[]> currentHeader = manifest.getHeaderContents();
        String[] lastRowInHeader = currentHeader.get(currentHeader.size() - 1); // "Loci Count" which needs to be last to terminate the header...
        for (int i = 0; i < currentHeader.size() - 1; i++) {
            String[] rowValues = currentHeader.get(i);
            if (numColumns == -1) {
                numColumns = rowValues.length;
            }
            addHeaderLine(output, numColumns, rowValues);
        }
        addHeaderLine(output, numColumns, Build37ExtendedIlluminaManifest.EXTENDED_MANIFEST_VERSION_HEADER_NAME, VERSION);
        addHeaderLine(output, numColumns, Build37ExtendedIlluminaManifest.EXTENDED_MANIFEST_TARGET_BUILD_HEADER_NAME, TARGET_BUILD);
        addHeaderLine(output, numColumns, Build37ExtendedIlluminaManifest.EXTENDED_MANIFEST_TARGET_REFERENCE_HEADER_NAME, TARGET_REFERENCE_FILE.getAbsolutePath());
        addHeaderLine(output, numColumns, Build37ExtendedIlluminaManifest.EXTENDED_MANIFEST_CLUSTER_FILE_HEADER_NAME, CLUSTER_FILE.getAbsolutePath());
        if (DBSNP_FILE != null) {
            addHeaderLine(output, numColumns, Build37ExtendedIlluminaManifest.EXTENDED_MANIFEST_DBSNP_FILE_HEADER_NAME, DBSNP_FILE.getAbsolutePath());
        }

        final String[] supportedBuildsFields = new String[SUPPORTED_BUILD.size() + 1];
        final String[] supportedReferenceFileFields = new String[SUPPORTED_BUILD.size() + 1];
        final String[] supportedChainFileFields = new String[SUPPORTED_BUILD.size() + 1];
        supportedBuildsFields[0] = Build37ExtendedIlluminaManifest.EXTENDED_MANIFEST_SUPPORTED_BUILD_HEADER_NAME;
        supportedReferenceFileFields[0] = Build37ExtendedIlluminaManifest.EXTENDED_MANIFEST_SUPPORTED_REFERENCE_HEADER_NAME;
        supportedChainFileFields[0] = Build37ExtendedIlluminaManifest.EXTENDED_MANIFEST_SUPPORTED_CHAIN_FILE_HEADER_NAME;
        for (int i = 0; i < SUPPORTED_BUILD.size(); i++) {
            supportedBuildsFields[i + 1] = SUPPORTED_BUILD.get(i);
            supportedReferenceFileFields[i + 1] = SUPPORTED_REFERENCE_FILE.get(i).getAbsolutePath();
            supportedChainFileFields[i + 1] = SUPPORTED_CHAIN_FILE.get(i).getAbsolutePath();
        }
        addHeaderLine(output, numColumns, supportedBuildsFields);
        addHeaderLine(output, numColumns, supportedReferenceFileFields);
        addHeaderLine(output, numColumns, supportedChainFileFields);
        addHeaderLine(output, numColumns, lastRowInHeader);

        addHeaderLine(output, numColumns, "[Assay]");

        // write the extended headers
        final String[] extendedHeader = ArrayUtils.addAll(manifest.getManifestFileHeaderNames(), Build37ExtendedIlluminaManifest.EXTENDED_MANIFEST_HEADERS);
        output.write(StringUtils.join(extendedHeader, ","));
        output.newLine();
    }

    private void addHeaderLine(final BufferedWriter out, final int numColumns, final String... fields) throws IOException {
        String[] rowValues = new String[numColumns];
        System.arraycopy(fields, 0, rowValues, 0, fields.length);
        out.write(StringUtils.join(rowValues, ","));
        out.newLine();
    }
}


