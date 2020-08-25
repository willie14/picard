package picard.arrays.illumina;

import picard.PicardException;

/**
 * A simple class to unify the mostly identical fields in Illumina's binary Bead Pool Manifest file (.bpm) and
 * it's text version.
 * This class also contains additional fields used by GtcToVcf to create a VCF from Illumina GTC files.
 */
public class ExtendedIlluminaManifestRecord {
    private final IlluminaBPMLocusEntry locusEntry;
    private final IlluminaManifestRecord illuminaManifestRecord;

    private Build37ExtendedIlluminaManifestRecord.Flag flag = Build37ExtendedIlluminaManifestRecord.Flag.PASS;




    public ExtendedIlluminaManifestRecord(final IlluminaBPMLocusEntry locusEntry,
                                          final IlluminaManifestRecord illuminaManifestRecord,
                                          final boolean dupe) {
        this.locusEntry = locusEntry;
        this.illuminaManifestRecord = illuminaManifestRecord;

        validate();

        //set dupe first so it can be overridden by fail flags
        if (dupe) flag = Build37ExtendedIlluminaManifestRecord.Flag.DUPE;

        // Look for entries which Illumina has marked as invalid
        if (getChrom().equals(IlluminaManifestRecord.ILLUMINA_FLAGGED_BAD_CHR)) {
            flag = Build37ExtendedIlluminaManifestRecord.Flag.ILLUMINA_FLAGGED;
        }
//
//        if (!r.getMajorGenomeBuild().trim().equals(BUILD_36) && !r.getMajorGenomeBuild().trim().equals(BUILD_37)) {
//            flag = Build37ExtendedIlluminaManifestRecord.Flag.UNSUPPORTED_GENOME_BUILD;
//        }


    }

    public Build37ExtendedIlluminaManifestRecord.Flag getFlag() {
        return flag;
    }

    public Boolean isDupe() {
        return flag == Build37ExtendedIlluminaManifestRecord.Flag.DUPE;
    }



    private void validate() {
        validateEntryField(locusEntry.ilmnId, illuminaManifestRecord.getIlmnId(), "ilmnId");
        validateEntryField(locusEntry.name, illuminaManifestRecord.getName(), "name");
        validateEntryField(locusEntry.ilmnStrand, illuminaManifestRecord.getIlmnStrand(), "ilmnStrand");
        validateEntryField(locusEntry.snp, illuminaManifestRecord.getSnp(), "snp");
        validateEntryField(locusEntry.chrom, illuminaManifestRecord.getChr(), "chrom");
        validateEntryField(locusEntry.ploidy, illuminaManifestRecord.getPloidy(), "ploidy");
        validateEntryField(locusEntry.species, illuminaManifestRecord.getSpecies(), "species");
        validateEntryField(locusEntry.mapInfo, illuminaManifestRecord.getPosition(), "mapInfo");
        validateEntryField(locusEntry.addressA, Integer.parseInt(illuminaManifestRecord.getAddressAId()), "addressAId");
        if (locusEntry.version == 4) {
            validateEntryField(locusEntry.alleleAProbeSeq, illuminaManifestRecord.getAlleleAProbeSeq(), "alleleAProbeSeq");
        }
        if ((locusEntry.addressB != -1) && (illuminaManifestRecord.getAddressBId() != null)) {
            validateEntryField(locusEntry.addressB, Integer.parseInt(illuminaManifestRecord.getAddressBId()), "addressBId");
        }
        if (locusEntry.version == 4) {
            validateEntryField(locusEntry.alleleBProbeSeq, illuminaManifestRecord.getAlleleBProbeSeq(), "alleleBProbeSeq");
        }
        validateEntryField(locusEntry.genomeBuild, illuminaManifestRecord.getGenomeBuild(), "genomeBuild");
        validateEntryField(locusEntry.source, illuminaManifestRecord.getSource(), "source");
        validateEntryField(locusEntry.sourceVersion, illuminaManifestRecord.getSourceVersion(), "sourceVersion");
        validateEntryField(locusEntry.sourceStrand, illuminaManifestRecord.getSourceStrand(), "sourceStrand");
        if (locusEntry.version == 4) {
            validateEntryField(locusEntry.sourceSeq, illuminaManifestRecord.getSourceSeq(), "sourceSeq");
            validateEntryField(locusEntry.topGenomicSeq, illuminaManifestRecord.getTopGenomicSeq(), "topGenomicSeq");
        }
//        if (record.getExpClusters() != null) {
//            validateEntryField(locusEntry.expClusters, Integer.parseInt(record.getExpClusters()), "expClusters");
//        }
//        validateEntryField(locusEntry.intensityOnly, record.getIntensityOnly(), "intensityOnly");
        if (locusEntry.version == 8) {
            validateEntryField(locusEntry.refStrand, illuminaManifestRecord.getRefStrand(), "refStrand");
        }
    }

    private String getChrom() {
        // TODO - previously I did a '.trim()
        return locusEntry.chrom;
    }

    private void validateEntryField(final Object locusEntryField, final Object recordField, final String fieldName) {
        // TODO - should this just be a non-passing entry?  Probably so.
        if (!locusEntryField.equals(recordField)) {
            throw new PicardException("Field '" + fieldName + "' disagrees between BPM file (found '" + locusEntryField + "') and CSV (found: '" + recordField + "')");
        }
    }
}
