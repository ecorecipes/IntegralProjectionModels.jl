"""
R expression lexer and recursive-descent parser.

Tokenizes R expression strings from PADRINO and parses them into an AST
(`RExpr`) that can be translated to Julia expressions.
"""

# --- Token types ---

@enum TokenType begin
    TOK_NUMBER
    TOK_IDENT
    TOK_STRING
    TOK_PLUS
    TOK_MINUS
    TOK_STAR
    TOK_SLASH
    TOK_CARET
    TOK_LPAREN
    TOK_RPAREN
    TOK_COMMA
    TOK_EQUALS
    TOK_LT
    TOK_GT
    TOK_LE
    TOK_GE
    TOK_EQ
    TOK_NE
    TOK_AND
    TOK_OR
    TOK_NOT
    TOK_COLON
    TOK_EOF
end

struct Token
    type::TokenType
    value::String
    pos::Int
end

# --- Lexer ---

function tokenize(expr::AbstractString)
    tokens = Token[]
    i = 1
    n = length(expr)

    while i <= n
        c = expr[i]

        # Skip whitespace
        if isspace(c)
            i += 1
            continue
        end

        # Numbers (including decimals and scientific notation)
        if isdigit(c) || (c == '.' && i + 1 <= n && isdigit(expr[i + 1]))
            start = i
            while i <= n && (isdigit(expr[i]) || expr[i] == '.')
                i += 1
            end
            # Scientific notation
            if i <= n && (expr[i] == 'e' || expr[i] == 'E')
                i += 1
                if i <= n && (expr[i] == '+' || expr[i] == '-')
                    i += 1
                end
                while i <= n && isdigit(expr[i])
                    i += 1
                end
            end
            push!(tokens, Token(TOK_NUMBER, expr[start:i-1], start))
            continue
        end

        # Identifiers (including dots for R names like is.na)
        if isletter(c) || c == '_' || c == '.'
            start = i
            while i <= n && (isletter(expr[i]) || isdigit(expr[i]) || expr[i] == '_' || expr[i] == '.')
                i += 1
            end
            push!(tokens, Token(TOK_IDENT, expr[start:i-1], start))
            continue
        end

        # String literals
        if c == '"' || c == '\''
            quote_char = c
            start = i
            i += 1
            while i <= n && expr[i] != quote_char
                if expr[i] == '\\' && i + 1 <= n
                    i += 1
                end
                i += 1
            end
            i += 1  # skip closing quote
            push!(tokens, Token(TOK_STRING, expr[start+1:i-2], start))
            continue
        end

        # Two-character operators
        if i + 1 <= n
            two = expr[i:i+1]
            if two == "<="
                push!(tokens, Token(TOK_LE, "<=", i)); i += 2; continue
            elseif two == ">="
                push!(tokens, Token(TOK_GE, ">=", i)); i += 2; continue
            elseif two == "=="
                push!(tokens, Token(TOK_EQ, "==", i)); i += 2; continue
            elseif two == "!="
                push!(tokens, Token(TOK_NE, "!=", i)); i += 2; continue
            elseif two == "&&"
                push!(tokens, Token(TOK_AND, "&&", i)); i += 2; continue
            elseif two == "||"
                push!(tokens, Token(TOK_OR, "||", i)); i += 2; continue
            end
        end

        # Single-character operators
        if c == '+'
            push!(tokens, Token(TOK_PLUS, "+", i))
        elseif c == '-'
            push!(tokens, Token(TOK_MINUS, "-", i))
        elseif c == '*'
            push!(tokens, Token(TOK_STAR, "*", i))
        elseif c == '/'
            push!(tokens, Token(TOK_SLASH, "/", i))
        elseif c == '^'
            push!(tokens, Token(TOK_CARET, "^", i))
        elseif c == '('
            push!(tokens, Token(TOK_LPAREN, "(", i))
        elseif c == ')'
            push!(tokens, Token(TOK_RPAREN, ")", i))
        elseif c == ','
            push!(tokens, Token(TOK_COMMA, ",", i))
        elseif c == '='
            push!(tokens, Token(TOK_EQUALS, "=", i))
        elseif c == '<'
            push!(tokens, Token(TOK_LT, "<", i))
        elseif c == '>'
            push!(tokens, Token(TOK_GT, ">", i))
        elseif c == '!'
            push!(tokens, Token(TOK_NOT, "!", i))
        elseif c == ':'
            push!(tokens, Token(TOK_COLON, ":", i))
        elseif c == '&'
            push!(tokens, Token(TOK_AND, "&", i))
        elseif c == '|'
            push!(tokens, Token(TOK_OR, "|", i))
        else
            error("Unexpected character '$c' at position $i in expression: $expr")
        end
        i += 1
    end

    push!(tokens, Token(TOK_EOF, "", n + 1))
    return tokens
end

# --- AST nodes ---

abstract type RExpr end

struct RNumber <: RExpr
    value::Float64
end

struct RIdent <: RExpr
    name::String
end

struct RString <: RExpr
    value::String
end

struct RCall <: RExpr
    func::String
    args::Vector{RExpr}
    kwargs::Vector{Pair{String, RExpr}}
end

struct RBinOp <: RExpr
    op::String
    left::RExpr
    right::RExpr
end

struct RUnaryOp <: RExpr
    op::String
    operand::RExpr
end

# --- Parser ---

mutable struct Parser
    tokens::Vector{Token}
    pos::Int
end

Parser(tokens::Vector{Token}) = Parser(tokens, 1)

function peek(p::Parser)
    return p.tokens[p.pos]
end

function advance!(p::Parser)
    t = p.tokens[p.pos]
    p.pos += 1
    return t
end

function expect!(p::Parser, type::TokenType)
    t = advance!(p)
    if t.type != type
        error("Expected $type but got $(t.type) ('$(t.value)') at position $(t.pos)")
    end
    return t
end

function parse_expr(tokens::Vector{Token})
    p = Parser(tokens)
    result = parse_or(p)
    return result
end

function parse_rexpr(expr::AbstractString)
    tokens = tokenize(expr)
    return parse_expr(tokens)
end

# Precedence levels (lowest to highest):
# or (||, |)
# and (&&, &)
# not (!)
# comparison (<, >, <=, >=, ==, !=)
# addition (+, -)
# multiplication (*, /)
# unary (-, +)
# power (^)
# function call / primary

function parse_or(p::Parser)
    left = parse_and(p)
    while peek(p).type == TOK_OR
        advance!(p)
        right = parse_and(p)
        left = RBinOp("||", left, right)
    end
    return left
end

function parse_and(p::Parser)
    left = parse_comparison(p)
    while peek(p).type == TOK_AND
        advance!(p)
        right = parse_comparison(p)
        left = RBinOp("&&", left, right)
    end
    return left
end

function parse_comparison(p::Parser)
    left = parse_addition(p)
    while peek(p).type in (TOK_LT, TOK_GT, TOK_LE, TOK_GE, TOK_EQ, TOK_NE)
        op = advance!(p).value
        right = parse_addition(p)
        left = RBinOp(op, left, right)
    end
    return left
end

function parse_addition(p::Parser)
    left = parse_multiplication(p)
    while peek(p).type in (TOK_PLUS, TOK_MINUS)
        op = advance!(p).value
        right = parse_multiplication(p)
        left = RBinOp(op, left, right)
    end
    return left
end

function parse_multiplication(p::Parser)
    left = parse_unary(p)
    while peek(p).type in (TOK_STAR, TOK_SLASH)
        op = advance!(p).value
        right = parse_unary(p)
        left = RBinOp(op, left, right)
    end
    return left
end

function parse_unary(p::Parser)
    if peek(p).type == TOK_MINUS
        advance!(p)
        operand = parse_power(p)
        return RUnaryOp("-", operand)
    elseif peek(p).type == TOK_NOT
        advance!(p)
        operand = parse_power(p)
        return RUnaryOp("!", operand)
    end
    return parse_power(p)
end

function parse_power(p::Parser)
    base = parse_primary(p)
    if peek(p).type == TOK_CARET
        advance!(p)
        # Right-associative
        exp = parse_unary(p)
        return RBinOp("^", base, exp)
    end
    return base
end

function parse_primary(p::Parser)
    t = peek(p)

    if t.type == TOK_NUMBER
        advance!(p)
        return RNumber(parse(Float64, t.value))
    end

    if t.type == TOK_STRING
        advance!(p)
        return RString(t.value)
    end

    if t.type == TOK_IDENT
        advance!(p)
        name = t.value

        # Check for function call
        if peek(p).type == TOK_LPAREN
            advance!(p)  # consume '('
            args = RExpr[]
            kwargs = Pair{String, RExpr}[]

            if peek(p).type != TOK_RPAREN
                _parse_call_arg!(p, args, kwargs)
                while peek(p).type == TOK_COMMA
                    advance!(p)  # consume ','
                    _parse_call_arg!(p, args, kwargs)
                end
            end
            expect!(p, TOK_RPAREN)
            return RCall(name, args, kwargs)
        end

        # Check for range operator ':'
        if peek(p).type == TOK_COLON
            advance!(p)
            right = parse_primary(p)
            return RCall("seq", [RIdent(name), right], Pair{String, RExpr}[])
        end

        return RIdent(name)
    end

    if t.type == TOK_LPAREN
        advance!(p)
        expr = parse_or(p)
        expect!(p, TOK_RPAREN)
        return expr
    end

    error("Unexpected token $(t.type) ('$(t.value)') at position $(t.pos)")
end

function _parse_call_arg!(p::Parser, args, kwargs)
    # Look ahead: is this a keyword argument? (name = expr)
    if peek(p).type == TOK_IDENT
        saved_pos = p.pos
        name_tok = advance!(p)
        if peek(p).type == TOK_EQUALS
            advance!(p)  # consume '='
            val = parse_or(p)
            push!(kwargs, name_tok.value => val)
            return
        end
        # Not a kwarg, rewind
        p.pos = saved_pos
    end
    push!(args, parse_or(p))
end
