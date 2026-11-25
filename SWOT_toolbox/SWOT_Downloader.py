"""Download SWOT Raster data for a bounding box and time interval."""
import argparse


def main():
    import earthaccess
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("destination")
    parser.add_argument("start", help="YYYY-MM-DD")
    parser.add_argument("end", help="YYYY-MM-DD")
    parser.add_argument("--bbox", nargs=4, type=float, required=True,
                        metavar=("WEST", "SOUTH", "EAST", "NORTH"))
    parser.add_argument("--collection", default="SWOT_L2_HR_Raster_D")
    args = parser.parse_args()
    earthaccess.login()
    results = earthaccess.search_data(short_name=args.collection,
                                      bounding_box=tuple(args.bbox),
                                      temporal=(args.start, args.end))
    earthaccess.download(results, args.destination)


if __name__ == "__main__":
    main()
