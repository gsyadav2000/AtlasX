from atlasx.io import BEDReader


def main():

    reader = BEDReader("data/example/example.bed")

    peaks = reader.load()

    print("=" * 50)
    print("BED Reader Test")
    print("=" * 50)

    print(f"Total Peaks : {len(peaks)}")
    print()

    print("First Peak:")
    print(peaks[0])

    print()

    peaks[0].summary()


if __name__ == "__main__":
    main()