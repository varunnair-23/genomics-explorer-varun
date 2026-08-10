metrics <- read.csv("results/dna_metrics.csv")

print(metrics)

g_count <- metrics$Count[metrics$Nucleotide == "G"]
c_count <- metrics$Count[metrics$Nucleotide == "C"]
total_count <- sum(metrics$Count)

gc_content <- ((g_count + c_count) / total_count) * 100

print(paste("GC-content percentage:", gc_content))

png(
  filename = "results/week1_chart.png",
  width = 800,
  height = 600
)

barplot(
  metrics$Count,
  names.arg = metrics$Nucleotide,
  main = paste("DNA Nucleotide Counts | GC-content:", round(gc_content, 2), "%"),
  xlab = "Nucleotide",
  ylab = "Count",
  col = c("gray20", "gray40", "gray60", "gray80")
)

dev.off()

print("Chart saved to results/week1_chart.png")