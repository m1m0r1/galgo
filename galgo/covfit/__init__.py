

def load_gender_info(gender_info):
    sample_genders = {}
    with open(gender_info) as fp:
        for line in fp:
            sample, gender = line.rstrip().split()
            gender = int(gender)
            assert gender in (1, 2)
            sample_genders[sample] = gender

    return sample_genders

def get_ploidies(chrom, genders):
    if chrom in ('X', 'chrX'):
        return tuple(1 if s == 1 else 2 for s in genders)
    elif chrom in ('Y', 'chrY'):
        return tuple(1 if s == 1 else 0 for s in genders)
    return (2,) * len(genders)
