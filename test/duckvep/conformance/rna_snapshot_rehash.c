#include "../../../src/core/duckvep_core_snapshot.c"

int main(int argc, char **argv)
{
	FILE *file;
	unsigned char *data;
	long size;
	uint64_t hash;
	snap_header_t *header;
	size_t section;

	if ((argc != 2 && argc != 3) || (file = fopen(argv[1], "rb")) == NULL)
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
	if (argc == 3) {
		char *end;
		uint64_t padded;
		section = strtoul(argv[2], &end, 10);
		if (*end != '\0' || section >= SNAP_SECTIONS ||
		    header->offset[section] > (uint64_t)size ||
		    header->bytes[section] > (uint64_t)size - header->offset[section])
			return 5;
		padded = (header->bytes[section] + SNAP_ALIGN - 1u) / SNAP_ALIGN * SNAP_ALIGN;
		if (padded > (uint64_t)size - header->offset[section])
			return 5;
		memmove(data + header->offset[section], data + header->offset[section] + padded,
		    (size_t)size - (size_t)header->offset[section] - (size_t)padded);
		if (section < SNAP_ARRAYS) header->present &= ~(UINT64_C(1) << section);
		header->bytes[section] = 0;
		while (++section < SNAP_SECTIONS)
			header->offset[section] -= padded;
		size -= (long)padded;
		header->file_bytes = (uint64_t)size;
	}
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
	if (fclose(file) != 0 || (file = fopen(argv[1], "wb")) == NULL ||
	    fwrite(data, 1, (size_t)size, file) != (size_t)size)
		return 6;
	free(data);
	return fclose(file) != 0;
}
