test_correlationFinder <- function()
{
	results <- correlationFinder()
	checkTrue(length(results$YOR063W) == 1)
	checkTrue(results$YOR063W == "YOR167C")
}