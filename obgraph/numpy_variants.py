import numpy as np
import logging


class NumpyVariants:
    properties = {"header", "variants"}
    def __init__(self, header, variants):
        self.header = header
        self.variants = variants

    def to_file(self, file_name):
        np.savez(file_name,
            header=self.header,
            variants=self.variants
        )

    @classmethod
    def from_file(cls, file_name):
        data = np.load(file_name)
        return cls(data["header"], data["variants"])

    def to_vcf_with_genotypes(self, file_name, sample_name, genotypes, add_header_lines=None, ignore_homo_ref=False, add_genotype_likelyhoods=None):
        logging.info("Writing to file %s" % file_name)
        if ignore_homo_ref:
            logging.info("Will not write variants with genotype 0/0")
        with open(file_name, "w") as f:
            # f.write(self._header_lines)
            for header_line in self.header:  # last element is empty
                if header_line.startswith("#CHROM"):
                    if sample_name != "":
                        header_line = header_line.strip() + "\t" + sample_name + "\n"
                    if add_header_lines is not None:
                        for additional_header_line in add_header_lines:
                            f.writelines([additional_header_line + "\n"])
                f.writelines([header_line])

            lines = []
            for i, (variant, genotype) in enumerate(zip(self.variants, genotypes)):
                genotype = genotype.decode("utf-8")
                if ignore_homo_ref and genotype == "0/0":
                    continue

                variant = variant.decode("utf-8").strip()
                if i % 1000000 == 0:
                    logging.info("%d variants written to file." % i)

                line = "%s\t%s\n" % (variant, genotype)
                if add_genotype_likelyhoods is not None:
                    likelyhoods = add_genotype_likelyhoods[i]
                    phred_likelyhoods = [
                        #float(round(np.maximum(0, np.minimum(255, -prob * np.log10(np.e))), 4))
                        int(np.minimum(2550000000000, -10 * prob * np.log10(np.e)))
                        for prob in likelyhoods
                    ]

                    genotype_likelyhoods = [
                        round(prob * np.log10(np.e), 0) for prob in likelyhoods
                    ]

                    phred_likelyhoods_str = ",".join(str(p) for p in phred_likelyhoods)
                    gl_str = ",".join(str(p) for p in genotype_likelyhoods)
                    #print(likelyhoods, phred_likelyhoods, np.exp(likelyhoods))
                    line = "%s:PL:GL\t%s:%s:%s\n" % (variant, genotype, phred_likelyhoods_str, gl_str)
                    #print(line)

                lines.append(line)

            f.writelines(lines)

    @classmethod
    def from_vcf(cls, vcf_file_name):
        header = []
        variants = []
        with open(vcf_file_name, "r") as f:
            for i, line in enumerate(f):
                #line_decoded = line.decode("utf-8")
                if line.startswith("#"):
                    header.append(line)
                else:
                    # add byte line
                    # remove info field to save space
                    l = line.split()
                    l[7] = "."
                    line = "\t".join(l)
                    variants.append(str.encode(line))  # add as byte

                if i % 10000 == 0:
                    logging.info("%d variants processed" % i)

        variants = np.array(variants)
        logging.info("Variants type: %s" % variants.dtype)
        return cls(np.array(header, dtype=str), variants)


