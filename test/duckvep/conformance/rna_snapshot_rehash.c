#include "../../../src/core/duckvep_core_snapshot.c"

int main(int argc, char **argv)
{
	FILE *file;
	unsigned char *data;
	long size;
	uint64_t hash;
	snap_header_t *header;
	size_t section;

	if (argc != 2 || (file = fopen(argv[1], "rb+")) == NULL)
		return 1;
	if (fseek(file, 0, SEEK_END) || (size = ftell(file)) < 0 ||
	    fseek(file, 0, SEEK_SET))
		return 2;
	data = malloc((size_t)size);
	if (data == NULL || fread(data, 1, (size_t)size, file) != (size_t)size)
		return 3;
	header = (snap_header_t *)data;
	if (header->version != SNAP_VERSION ||
	    header->header_bytes != sizeof(*header))
		return 4;
	header->checksum = 0;
	hash = snap_hash(header->version, data, sizeof(*header));
	for (section = 0; section < SNAP_SECTIONS; section++) {
		if (header->offset[section] > (uint64_t)size ||
		    header->bytes[section] > (uint64_t)size - header->offset[section])
			return 5;
		hash = snap_hash(hash, data + header->offset[section],
		    (size_t)header->bytes[section]);
	}
	header->checksum = hash;
	if (fseek(file, 0, SEEK_SET) ||
	    fwrite(data, 1, (size_t)size, file) != (size_t)size)
		return 6;
	free(data);
	return fclose(file) != 0;
}
