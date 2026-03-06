# LLM Provider System - Single file for all LLM functionality

#' Call LLM API
#'
#' Routes prompts to OpenAI, Anthropic, Google Gemini, or local Ollama models via ellmer package.
#'
#' @param prompt Character. The prompt to send.
#' @param provider Character. One of "openai", "anthropic", "gemini", "ollama", or "vllm".
#' @param model Character. Model name (e.g., "gpt-5.2", "claude-opus-4-6", "gemini-2.0-flash").
#'   For Ollama, use model names like "llama2", "mistral", "neural-chat", etc.
#' @param params An ellmer params object created by ellmer::params().
#' @param api_key Character. API key (uses env var if NULL).
#'   OpenAI: OPENAI_API_KEY
#'   Anthropic: ANTHROPIC_API_KEY
#'   Gemini: GOOGLE_API_KEY or GEMINI_API_KEY
#'   Ollama/vLLM: Not required (local)
#' @param system_prompt Character. System prompt.
#' @param api_url Character. API URL for local LLMs (Ollama or vLLM).
#'   Ollama default: "http://localhost:11434"
#'   vLLM default: "http://localhost:8000"
#'
#' @details
#' **Ollama Setup:**
#' To use Ollama locally, first install it from https://ollama.ai, then run:
#' ```
#' ollama pull llama2       # or another model
#' ollama serve             # start the server (runs on localhost:11434)
#' ```
#' Popular models: llama2, mistral, neural-chat, dolphin-mixtral, etc.
#'
#' **vLLM Setup:**
#' Similar to Ollama but optimized for performance. Requires Python installation.
#'
#' @return Character. The LLM response text.
#' @importFrom ellmer params chat_openai chat_anthropic chat_google_gemini chat_ollama chat_vllm
#' @export
call_llm <- function(prompt, provider = "openai", model = NULL,
                     params = NULL, api_key = NULL, system_prompt = NULL,
                     api_url = NULL) {
  # Ollama models often include a tag suffix (e.g. "llama3.2:latest").
  # The ellmer client expects the bare model name, so strip any tag.
  if (!is.null(model) && provider == "ollama") {
    model <- sub(":.*$", "", model)
  }
  if (!requireNamespace("ellmer", quietly = TRUE)) {
    stop("Package 'ellmer' required but not installed")
  }
  # warn if provider uses API key and none provided
  if (provider %in% c("openai", "anthropic", "gemini")) {
    envs <- switch(provider,
                   openai = "OPENAI_API_KEY",
                   anthropic = "ANTHROPIC_API_KEY",
                   gemini = c("GOOGLE_API_KEY", "GEMINI_API_KEY"))
    if (is.null(api_key)) {
      has <- vapply(envs, function(e) nzchar(Sys.getenv(e, "")), logical(1))
      if (!any(has)) {
        message("[GPTAnno] warning: no API key found for provider '", provider, "'. ",
                "Set one via environment variable (", paste(envs, collapse = "/"), ") ",
                "or pass via llm_config$api_key.")
      }
    }
  }

  if (provider == "openai") {
    chat <- ellmer::chat_openai(
      system_prompt = system_prompt,
      model = model,
      params = params,
      api_key = api_key,
      echo = "none"
    )
  } else if (provider == "anthropic") {
    chat <- ellmer::chat_anthropic(
      system_prompt = system_prompt,
      model = model,
      params = params,
      api_key = api_key,
      echo = "none"
    )
  } else if (provider == "gemini") {
    chat <- ellmer::chat_google_gemini(
      system_prompt = system_prompt,
      model = model,
      params = params,
      api_key = api_key,
      echo = "none"
    )
  } else if (provider == "ollama") {
    # Set default Ollama API URL if not provided
    if (is.null(api_url)) {
      api_url <- "http://localhost:11434"
    }
    chat <- ellmer::chat_ollama(
      system_prompt = system_prompt,
      model = model,
      params = params,
      base_url = api_url,
      echo = "none"
    )
  } else if (provider == "vllm") {
    # Set default vLLM API URL if not provided
    if (is.null(api_url)) {
      api_url <- "http://localhost:8000"
    }
    chat <- ellmer::chat_vllm(
      system_prompt = system_prompt,
      model = model,
      params = params,
      base_url = api_url,
      echo = "none"
    )
  } else {
    stop("Unknown provider: ", provider, ". Use 'openai', 'anthropic', 'gemini', 'ollama', or 'vllm'")
  }

  return(chat$chat(prompt))
}

#' Prepare LLM Configuration
#'
#' Helper function that merges legacy `model` parameter with new `llm_config` structure.
#' Provides backward compatibility while supporting new multi-provider configuration.
#'
#' @param model Character. Model name (default: "gpt-5.2"). Used if llm_config not provided.
#' @param llm_config Optional list. Configuration with: provider, model, temperature, max_tokens, api_key, api_url, system_prompt.
#'
#' @return A list with merged LLM configuration.
#'
#' @keywords internal
#' @export
prepare_config <- function(model = "gpt-5.2", llm_config = NULL) {
  config <- list(
    provider = "openai",
    model = model,
    temperature = NULL,
    max_tokens = NULL,
    api_key = NULL,
    api_url = NULL,
    system_prompt = NULL
  )

  if (!is.null(llm_config)) {
    for (key in names(llm_config)) {
      config[[key]] <- llm_config[[key]]
    }
  }

  return(config)
}

# Helper: NULL coalescing
`%||%` <- function(x, y) if (is.null(x)) y else x

#' Check Ollama Service Availability
#'
#' Checks if Ollama is running and accessible at the specified URL.
#'
#' @param api_url Character. Ollama API URL (default: "http://localhost:11434")
#'
#' @return Logical. TRUE if Ollama is accessible, FALSE otherwise.
#'
#' @details
#' This function attempts to ping the Ollama API to verify it is running.
#' Returns FALSE silently if the service is not available.
#'
#' @examples
#' \dontrun{
#' check_ollama_available()
#' }
#'
#' @export
check_ollama_available <- function(api_url = "http://localhost:11434") {
  tryCatch({
    response <- httr::GET(paste0(api_url, "/api/tags"), httr::timeout(2))
    httr::status_code(response) == 200
  }, error = function(e) {
    FALSE
  })
}

#' List Available Ollama Models
#'
#' Lists all models currently available in Ollama.
#'
#' @param api_url Character. Ollama API URL (default: "http://localhost:11434")
#'
#' @return Character vector of model names, or NULL if Ollama is not available.
#'
#' @details
#' Requires Ollama to be running. If Ollama is not accessible, returns NULL with a warning.
#'
#' @examples
#' \dontrun{
#' list_ollama_models()
#' # Returns something like:
#' # [1] "llama2" "mistral" "neural-chat"
#' }
#'
#' @export
list_ollama_models <- function(api_url = "http://localhost:11434") {
  if (!check_ollama_available(api_url)) {
    warning("Ollama is not running at ", api_url, ". ",
            "Start Ollama with: ollama serve")
    return(NULL)
  }
  
  tryCatch({
    response <- httr::GET(paste0(api_url, "/api/tags"))
    content <- httr::content(response, as = "parsed")
    if (!is.null(content$models)) {
      return(sapply(content$models, function(x) x$name))
    }
    return(NULL)
  }, error = function(e) {
    warning("Error fetching Ollama models: ", e$message)
    return(NULL)
  })
}

# Package initialization
.onLoad <- function(libname, pkgname) {
  # Check and warn if ellmer is not installed
  if (!requireNamespace("ellmer", quietly = TRUE)) {
    packageStartupMessage(
      "GPTAnno requires the 'ellmer' package for LLM support. ",
      "Install with: install.packages('ellmer')"
    )
  }
  invisible(NULL)
}
