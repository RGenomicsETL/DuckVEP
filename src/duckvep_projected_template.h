/* Small-variant projected facts for presentation from one annotation cursor. */
static const char *const duckvep_projected_parts[] = {
    "WITH parameters AS (SELECT CAST(__DUCKVEP_MODEL__ AS VARCHAR) AS model_name, ",
    "CAST(__DUCKVEP_UPSTREAM__ AS UBIGINT) AS upstream_distance, ",
    "CAST(__DUCKVEP_DOWNSTREAM__ AS UBIGINT) AS downstream_distance), ",
    "source AS (SELECT e.*, ",
    "duckvep_allele_geometry(e.position, e.reference, e.alternate) AS geometry ",
    "FROM __DUCKVEP_EVENTS__ e), ",
    "projected AS (SELECT e.*, unnest(_duckvep_annotate_small_projected( ",
    "p.model_name, e.seq_region, e.position, e.reference, ",
    "CASE WHEN e.event_index IS NULL OR e.alternate = '<*>' ",
    "THEN error('duckvep_annotate_projected: expected literal small-variant events') ",
    "ELSE e.alternate END, ",
    "p.upstream_distance, p.downstream_distance)) AS projection ",
    "FROM source e CROSS JOIN parameters p) ",
    "SELECT * EXCLUDE (projection), projection.* FROM projected"
};
