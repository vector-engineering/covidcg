
with start_count as (
	select
		range_start,
		"count",
		sum("count") over (order by range_start asc rows between unbounded preceding and current row)
	from (
		select range_start, count(*) 
		from dna_coverage
		where reference = 'KX858756.1'
		group by range_start
	) start_sum
	order by range_start asc
), end_count as (
	select
		range_end,
		"count",
		sum("count") over (order by range_end asc rows between unbounded preceding and current row)
	from (
		select range_end, count(*) 
		from dna_coverage
		where reference = 'KX858756.1'
		group by range_end
	) end_sum
	order by range_end asc
)
select 
	start_count.range_start,
	end_count.range_end,
	start_count."sum" as "start_sum",
	end_count."sum" as "end_sum"
from start_count
full join end_count on start_count.range_start = end_count.range_end



with selected as (
	select sequence_id
	from metadata
	where region = 0
), start_count as (
	select range_start, count(*) 
	from dna_coverage
	inner join selected on dna_coverage.sequence_id = selected.sequence_id
	where reference = 'KX858756.1'
	group by range_start
), end_count as (
	select range_end, count(*) 
	from dna_coverage
	inner join selected on dna_coverage.sequence_id = selected.sequence_id
	where reference = 'KX858756.1'
	group by range_end
), start_end_count as (
	select
		coalesce(start_count.range_start, end_count.range_end) as "ind",
		coalesce(start_count."count", 0) as "start",
		coalesce(end_count."count", 0) as "end"
	from start_count
	full join end_count on start_count.range_start = end_count.range_end
	order by ind asc
)
select
	"ind",
	sum("start") over (order by "ind" asc) -
	sum("end") over (order by "ind" asc) as "diff"
from start_end_count
order by "ind" asc

with selected as (
	select sequence_id
	from metadata
	where region = 0
), start_count as (
	select gene, range_start, count(*)
	from gene_aa_coverage
	inner join selected on gene_aa_coverage.sequence_id = selected.sequence_id
	where reference = 'NC_001781.1' and gene in ('F', 'G')
	group by gene, range_start
), end_count as (
	select gene, range_end, count(*)
	from gene_aa_coverage
	inner join selected on gene_aa_coverage.sequence_id = selected.sequence_id
	where reference = 'NC_001781.1' and gene in ('F', 'G')
	group by gene, range_end
), start_end_count as (
	select
		coalesce(start_count.gene, end_count.gene) as "gene",
		coalesce(start_count.range_start, end_count.range_end) as "ind",
		coalesce(start_count."count", 0) as "start",
		coalesce(end_count."count", 0) as "end"
	from start_count
	full join end_count on 
		start_count.gene = end_count.gene and 
		start_count.range_start = end_count.range_end
	order by gene, ind asc
)
select
	"gene", "ind",
	sum("start") over (partition by "gene" order by "ind" asc) -
	sum("end") over (partition by "gene" order by "ind" asc) as "diff"
from start_end_count
order by "gene", "ind" asc